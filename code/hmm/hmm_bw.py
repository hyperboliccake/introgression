import numpy as np
from typing import List, Dict, Tuple


class HMM:
    def __init__(self):
        self.hidden_states = []
        self.observed_states = []
        self.hidden_states = []
        self.transitions = []
        self.emissions = []
        self.observations = []
        self.initial_p = []
        self.symbol_to_ind = {}

    def set_hidden_states(self, states: List[str]) -> None:
        '''
        Sets the hidden states of the HMM to the supplied list of strings
        '''
        self.hidden_states = states

    def set_observed_states(self, states: List[str]) -> None:
        '''
        Sets the observed states of the HMM to the supplied list of strings
        If not supplied will set to list of keys provided by emissions 
        '''
        self.observed_states = states

    def set_transitions(self, transitions: List[List[float]]) -> None:
        '''
        Set transition probabilities to the supplied list of values
        Element i,j is the probability of transitioning from state i to j
        Checks that each row sums to approximately 1
        '''
        self.transitions = np.array(transitions)
        # find rows which are not close to 1 as boolean array
        row_sum = np.sum(self.transitions, axis=1)
        out_of_range = np.invert(np.isclose(row_sum, 1))
        if np.any(out_of_range):
            # get first index out of range
            index = np.where(out_of_range)[0][0]
            raise ValueError(f"{transitions[index]} {row_sum[index]}")

    def set_emissions(self, emissions: List[Dict[str, float]]) -> None:
        '''
        Set the emission matrix of this HMM to the supplied values
        Input is list of dictionaries, index matches hidden state, keys match
            the states.
        Must have states set prior to calling
        '''
        # assumes all dicts have same keys
        if self.observed_states == []:
            self.observed_states = list(emissions[0].keys())
        self.symbol_to_ind = {key: ind
                              for ind, key in enumerate(self.observed_states)}

        self.emissions = np.array([[d[key] if key in d else 0
                                    for key in self.observed_states]
                                   for d in emissions])
        row_sum = np.sum(self.emissions, axis=1)
        out_of_range = np.invert(np.isclose(row_sum, 1))
        if np.any(out_of_range):
            # get first index out of range
            index = np.where(out_of_range)[0][0]
            raise ValueError(f"{self.observed_states[index]} {row_sum[index]}")

    def set_observations(self, observations: List[List[str]]) -> None:
        '''
        Set observations of the HMM.
        Input is a list of of sequences of states
        Must have states set prior to calling
        '''
        # one sequence for Viterbi, or list of observed sequences
        # for Baum-Welch; TODO maybe change this weirdness?
        if self.symbol_to_ind == []:
            raise AttributeError("Observations can not be\
                                 set without emissions")

        self.observations = np.array([[self.symbol_to_ind[s] for s in sequence]
                                      for sequence in observations])

    def set_initial_p(self, initial_p: List[float]) -> None:
        '''
        Set initial probabilities of the hidden states
        Checks that probability sums to 1
        '''

        self.initial_p = np.array(initial_p)
        assert np.isclose(sum(initial_p), 1), f"{initial_p} {sum(initial_p)}"

    def print_results(self, iterations: int, LL: float) -> None:
        '''
        Write current state of HMM to stdout
        '''
        print(
            f'''Iterations: {iterations}

Log Likelihood:
{LL:.30e}

Initial State Probabilities:'''
        )
        for i in range(len(self.hidden_states)):
            print(f"{self.hidden_states[i]}={self.initial_p[i]:.30e}")
        print()
        print("Transition Probabilities:")
        for i in range(len(self.hidden_states)):
            for j in range(len(self.hidden_states)):
                print(f"{self.hidden_states[i]},{self.hidden_states[j]}\
                    ={self.transitions[i][j]:.30e}")
        print()
        print("Emission Probabilities:")
        for i in range(len(self.hidden_states)):
            for k in sorted(self.observed_states):
                print(f"{self.hidden_states[i]},{k}=\
                      {self.emissions[i, self.symbol_to_ind[k]]:.30e}")
        print()

    def train(self,
           improvement_frac: float = .01,
           max_iterations: int = None) -> None:
        '''
        Train the hmm until either the max iterations is reached or
        the log likelihood fails to improve beyond the improvement factor
        '''

        # calculate current log likelihood
        print("calculating alpha")
        alpha = self.forward()

        LL = self.log_likelihood(alpha)

        iterations = 0
        self.print_results(iterations, LL)

        # continue until log likelihood has stopped increasing much
        threshold = improvement_frac * abs(LL)
        # to execute first iteration
        prev_LL = np.NINF
        while (max_iterations is not None
               and iterations < max_iterations)\
                or LL - prev_LL > threshold:

            print(f"Iteration {iterations}")

            print("calculating beta")
            beta = self.backward()
            print("calculating gamma")
            gamma = self.state_probs(alpha, beta)
            print("calculating xi")
            xi = self.bw(alpha, beta)

            print("updating parameters")

            self.initial_p = self.initial_probabilities(gamma)
            self.transitions = self.transition_probabilities(xi, gamma)
            self.emissions = self.emission_probabilities(gamma)

            assert np.isclose(np.sum(self.initial_p), 1), \
                f"{beta}\n{np.sum(self.initial_p)} {self.initial_p}"
            for t in self.transitions:
                assert np.isclose(sum(t), 1), \
                    f"{xi} {gamma} {sum(t)} {t}"
            for e in self.emissions:
                assert np.isclose(sum(e), 1), f"{sum(e.values())} {e}"

            iterations += 1

            print("calculating alpha")
            alpha = self.forward()

            prev_LL = LL
            LL = self.log_likelihood(alpha)

            # print results for every iteration
            self.print_results(iterations, LL)

            if LL < prev_LL and not np.isclose(LL, prev_LL):
                # NOTE does not stop execution
                print('PROBLEM: log-likelihood stopped increasing; \
                      stopping training now')

        print(f"finished in {iterations} iterations")

    def log_likelihood(self, alpha: np.array) -> float:
        '''
        Determines the log likelihood from the supplied observation of alpha
        '''
        # for multiple observations, product of LL
        return np.sum(
            np.logaddexp.reduce(
                alpha[:, -1, :],
                axis=1)
        )

    def transition_probabilities(self,
                                 xi: np.array,
                                 gamma: np.array) -> np.array:
        '''
        Caclulate new transition probabilities from xi and gamma
        '''
        # reduce along observation string
        num = np.logaddexp.reduce(xi, axis=1)
        num = np.logaddexp.reduce(num, axis=0)  # reduce along observations
        # reduce along observations, remove last
        den = np.logaddexp.reduce(gamma[:, :-1, :], axis=1)
        den = np.logaddexp.reduce(den, axis=0)  # reduce along observations

        transitions = np.exp(num - den[:, None])

        return transitions

    def initial_probabilities(self, gamma: np.array) -> np.array:
        '''
        Calculates the probability of beginning in each hidden state
        '''
        return np.exp(  # convert log likelihood to p
            np.logaddexp.reduce(gamma[:, 0, :])  # sum along states
        ) / len(self.observations)

    def emission_probabilities(self, gamma: np.array) -> np.array:
        '''
        Calcualte new emission probability matrix from gamma
        '''
        # denominator is sum of gamma along sequences and sequence length
        denominator = np.logaddexp.reduce(  # along sequences
            np.logaddexp.reduce(  # along sequence length
                gamma))

        # numerator is also sum along gamma but only where position ==
        # the observed symbols
        # convert observation to boolean matrix where
        # numerator[i, j, k] is true iff observed symbol at [j, k] == i
        numerator = np.array([i == self.observations
                              for i in range(len(self.observed_states))])
        # move axis to match gamma shape
        numerator = np.moveaxis(numerator, [0, 1, 2], [2, 0, 1])
        # use where to select gamma where numerator is true, NINF where false
        # broadcast to create a hidden_states x observed_states matrix
        numerator = np.where(numerator[:, :, None, :],
                             gamma[:, :, :, None],  # if true
                             np.NINF)  # if false
        # reduce along first two axes
        numerator = np.logaddexp.reduce(  # along sequences
            np.logaddexp.reduce(  # along sequence length
                numerator))

        return np.exp(numerator - denominator[:, None])

    def bw(self, alpha: np.array, beta: np.array) -> np.array:
        '''
        Calculate Baum-Welch xi from alpha and beta
        '''

        # shift times to match appropriate elements
        alpha = alpha[:, :-1, :]
        beta = beta[:, 1:, :]

        # tranpose emissions to make it easier to index into, take logs
        emis = np.transpose(self.emissions)
        emis = np.log(emis[self.observations[:, 1:]])

        # calculate numerator with broadcasting to proper size
        numerator = alpha[:, :, :, None] + beta[:, :, None, :] +\
            np.log(self.transitions) + emis[:, :, None, :]

        # denominator is summed along states
        denominator = np.logaddexp.reduce(
            np.logaddexp.reduce(numerator, axis=2),
            axis=2)

        xi = numerator - denominator[:, :, None, None]
        return xi

    def state_probs(self, alpha: np.array, beta: np.array) -> np.array:
        '''
        Calculate gamma from supplied alpha and beta
        '''
        # probability of being at state i at time j

        alpha_beta = alpha + beta
        denominator = np.logaddexp.reduce(alpha_beta, axis=2)
        return alpha_beta - denominator[:, :, None]

    def forward(self) -> np.array:
        '''
        Calculate alpha
        '''

        # probability that the sequence from 0 to t was observed and
        # Markov process was at state j at time t
        # returns array of size observations, observations[0], hidden_states
        # determine emission probabilities for each measured value
        emis = np.transpose(np.log(self.emissions[:, self.observations]))
        trans = np.log(self.transitions)
        alpha = np.empty((len(self.observations),
                          len(self.observations[0]),
                          len(self.hidden_states)), float)

        # initialize to initial probabilitiy * observed emission
        alpha[:, 0, :] = np.log(self.initial_p[None, :]) + emis[0, :, :]
        # recursively fill array
        for i in range(1, len(self.observations[0])):
            alpha[:, i, :] = np.logaddexp.reduce(alpha[:, i-1, :][:, :, None] +
                                                 trans[None, :, :], axis=1) \
                + emis[i, :, :]
        return alpha

    def backward(self) -> np.array:
        '''
        Calculate beta
        '''
        # probability that the sequence from t+1 to end was observed
        # and Markov process was at state j at time t
        emis = np.transpose(np.log(self.emissions[:, self.observations]))
        trans = np.log(self.transitions)
        beta = np.zeros((len(self.observations),
                         len(self.observations[0]),
                         len(self.hidden_states)), float)

        for i in range(len(self.observations[0]) - 2, -1, -1):
            beta[:, i, :] = np.logaddexp.reduce(
                (beta[:, i+1, :] + emis[i+1, :, :])[:, None, :] +
                trans[None, :, :], axis=2)

        return beta

    def calculate_max_states(self) -> Tuple[np.array, np.array]:
        '''
        Find the maximum likelihood hidden states and the corresponding 
        log probability for each state.
        Returned tuple is (probability, states)
        '''
        probabilities = np.empty((len(self.observations),
                                  len(self.hidden_states)), float)
        states = np.empty((len(self.observations),
                           len(self.hidden_states)), int)

        # build array of emissions based on observations
        emissions = np.log(np.transpose(self.emissions)[self.observations])

        trans_emis = np.log(self.transitions[None, :, :]) +\
            emissions[:, None, :]

        probabilities[0, :] = np.log(self.initial_p) + emissions[0]
        states[0, :] = -1

        for i in range(1, len(emissions)):
            prob = probabilities[i-1, :] + np.transpose(trans_emis[i])
            max_pos = np.argmax(prob, axis=1)
            probabilities[i, :] = [prob[j, a] for j, a in enumerate(max_pos)]
            states[i, :] = max_pos

        return probabilities, states

    def max_path(self, probs: np.array, states: np.array) -> np.array:
        '''
        Return the most likeli states for the HMM
        '''

        path = np.empty(len(self.observations), int)
        # state with most likely state
        path[0] = np.argmax(probs[-1])
        # fill in with most likely sequence from back
        for i, state in enumerate(reversed(states[1:])):
            path[i+1] = state[path[i]]

        # reverse to proper order
        path = np.flip(path)

        return path

    def viterbi(self) -> np.array:
        '''
        Find the maximum likelihood path for the first observation
        '''
        # consider only the first observation
        self.observations = self.observations[0]
        probs, states = self.calculate_max_states()
        return self.max_path(probs, states)

    def posterior_decoding(self) -> np.array:
        '''
        Return gamma as probability values
        probability[i, j, k] is for observation i, position j in observation
        and hidden state k
        '''
        alpha = self.forward()
        beta = self.backward()
        return np.exp(self.state_probs(alpha, beta))
