import math
from collections import defaultdict
import numpy as np
from typing import List, Dict

LOGZERO = 'LOGZERO'


def elog(n):
    if n == 0 or n == LOGZERO:
        return LOGZERO
    return math.log(n)


def elogproduct(m, n):
    if m == LOGZERO or n == LOGZERO:
        return LOGZERO
    return m + n


def elogsum(m, n):
    # approx(np.logaddexp(m, n))
    if m == LOGZERO:
        return n
    if n == LOGZERO:
        return m
    if m > n:
        return m + elog(1 + eexp(n - m))
    return n + elog(1 + eexp(m - n))


def eexp(n):
    if n == LOGZERO:
        return 0
    return math.e ** n


def ecompare(m, n):
    if m == LOGZERO:
        return False
    if n == LOGZERO:
        return True
    return m > n


class HMM:
    def __init__(self):
        # TODO replace with numpy arrays as much as possible

        self.states = []
        self.transitions = []
        self.emissions = []
        self.observations = []
        self.initial_p = []
        self.symbol_to_ind = {}

    def set_states(self, states: List[str]) -> None:
        '''
        Sets the states of the HMM to the supplied list of strings
        '''
        self.states = states

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
        self.symbols = list(emissions[0].keys())
        self.symbol_to_ind = {key: ind for ind, key in enumerate(self.symbols)}

        self.emissions = np.array([[d[key] for key in self.symbols]
                                   for d in emissions])
        row_sum = np.sum(self.emissions, axis=1)
        out_of_range = np.invert(np.isclose(row_sum, 1))
        if np.any(out_of_range):
            # get first index out of range
            index = np.where(out_of_range)[0][0]
            raise ValueError(f"{self.symbols[index]} {row_sum[index]}")

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

    def print_results(self, num_its: int, LL: float) -> None:
        '''
        Display current state of HMM to stdout
        '''
        print(
            f'''Iterations: {num_its}

Log Likelihood:
{LL:.30e}

Initial State Probabilities:'''
        )
        for i in range(len(self.states)):
            print(f"{self.states[i]}={self.initial_p[i]:.30e}")
        print()
        print("Transition Probabilities:")
        for i in range(len(self.states)):
            for j in range(len(self.states)):
                print(f"{self.states[i]},{self.states[j]}\
                    ={self.transitions[i][j]:.30e}")
        print()
        print("Emission Probabilities:")
        for i in range(len(self.states)):
            for k in sorted(self.symbols):
                print(f"{self.states[i]},{k}=\
                      {self.emissions[i, self.symbol_to_ind[k]]:.30e}")
        print()

    def go(self, improvement_frac=.01, num_its=0, prev_LL=0):

        # calculate current log likelihood
        print("calculating alpha")
        alpha = self.forward()

        # for multiple observations, product of LL
        # (note not base 2 log)
        LL_all = []
        LL = math.log(1)  # because we're taking product
        for seq in range(len(self.observations)):
            LL_seq = LOGZERO  # because we're taking sum
            for i in range(len(self.states)):
                LL_seq = elogsum(LL_seq, alpha[seq][-1][i])
            LL_all.append(LL_seq)
            LL = elogproduct(LL, LL_seq)
        self.print_results(num_its, LL)

        # continue until log likelihood has stopped increasing much
        threshold = improvement_frac * abs(LL)
        while num_its < 1 or LL - prev_LL > threshold:

            print('Iteration',  num_its)

            print("calculating beta")
            beta = self.backward()
            print("calculating gamma")
            gamma = self.state_probs(alpha, beta)
            print("calculating xi")
            xi = self.bw(alpha, beta)

            print("updating parameters")

            # initial probabilities
            pi = [LOGZERO for i in range(len(self.states))]
            for i in range(len(self.states)):
                for seq in range(len(self.observations)):
                        pi[i] = elogsum(pi[i], gamma[seq][0][i])
            pi = [eexp(x) / len(self.observations) for x in pi]

            # transition probabilities
            a = []
            for i in range(len(self.states)):
                row = []
                for j in range(len(self.states)):
                        num = LOGZERO
                        den = LOGZERO
                        for seq in range(len(self.observations)):
                            num_seq = LOGZERO
                            den_seq = LOGZERO
                            for o in range(len(self.observations[seq]) - 1):
                                # xi is probability of probability
                                # of being at state i at time t and
                                # state j at time t+1
                                num_seq = elogsum(num_seq, xi[seq][o][i][j])
                                # gamma is probability of being in
                                # state i at time t
                                den_seq = elogsum(den_seq, gamma[seq][o][i])
                            # weight numerator and denominator for the
                            # current observation sequence by 1/P(seq | model)
                            # num_seq = elogproduct(num_seq, -LL_all[seq])
                            # den_seq = elogproduct(den_seq, -LL_all[seq])
                            # add the current sequence contribution to total
                            num = elogsum(num, num_seq)
                            den = elogsum(den, den_seq)
                        assert den != LOGZERO, \
                            'probably something wrong with \
                            initial parameter values'
                        row.append(eexp(elogproduct(num, -den)))
                a.append(row)

            # emission probabilities
            b = []
            for state in range(len(self.states)):
                d = defaultdict(float)
                for symbol in range(len(self.emissions[state])):
                        num = LOGZERO
                        den = LOGZERO
                        for seq in range(len(self.observations)):
                            num_seq = LOGZERO
                            den_seq = LOGZERO
                            for o in range(len(self.observations[seq])):
                                if self.observations[seq][o] == symbol:
                                    # gamma is probability of
                                    # being in state i at time t
                                    num_seq = elogsum(num_seq,
                                                      gamma[seq][o][state])
                                den_seq = elogsum(den_seq,
                                                  gamma[seq][o][state])
                            # weight numerator and denominator for the
                            # current observation sequence by 1/P(seq | model)
                            # num_seq = elogproduct(num_seq, -LL_all[seq])
                            # den_seq = elogproduct(den_seq, -LL_all[seq])
                            # add the current sequence contribution to total
                            num = elogsum(num, num_seq)
                            den = elogsum(den, den_seq)
                        assert den != LOGZERO, \
                            'probably something wrong with \
                            initial parameter values'
                        d[self.symbols[symbol]] = eexp(elogproduct(num, -den))
                b.append(d)

            self.initial_p = pi
            self.transitions = a
            # TODO replace with just keeping matrix when not converting to dict
            self.set_emissions(b)
            # self.emissions = b

            assert np.isclose(sum(self.initial_p), 1), \
                f"{sum(self.initial_p)} {self.initial_p}"
            for t in self.transitions:
                assert np.isclose(sum(t), 1), f"{sum(t)} {t}"
            for e in self.emissions:
                assert np.isclose(sum(e), 1), f"{sum(e.values())} {e}"

            num_its += 1

            print("calculating alpha")
            alpha = self.forward()

            prev_LL = LL

            LL_all = []
            LL = math.log(1)  # because we're taking product
            for seq in range(len(self.observations)):
                LL_seq = LOGZERO  # because we're taking sum
                for i in range(len(self.states)):
                        LL_seq = elogsum(LL_seq, alpha[seq][-1][i])
                LL_all.append(LL_seq)
                LL = elogproduct(LL, LL_seq)

            # print results for every iteration
            self.print_results(num_its, LL)

            if LL < prev_LL and not np.isclose(LL, prev_LL):
                # NOTE does not stop execution
                print('PROBLEM: log-likelihood stopped increasing; \
                      stopping training now')

        print(f"finished in {num_its} iterations")

    def bw(self, alpha, beta):

        # probability of being at state i at time t and state j at
        # time t+1
        xi = []
        for seq in range(len(self.observations)):
            xi_current = []
            for o in range(len(self.observations[seq]) - 1):
                norm = LOGZERO
                matrix = [[LOGZERO] * len(self.states)
                          for i in range(len(self.states))]
                for i in range(len(self.states)):
                        for j in range(len(self.states)):
                            prob = elogproduct(
                                elog(self.emissions[j][self.observations[seq][o+1]]),
                                beta[seq][o+1][j])
                            prob = elogproduct(elog(self.transitions[i][j]), prob)
                            prob = elogproduct(alpha[seq][o][i], prob)
                            matrix[i][j] = prob
                            norm = elogsum(norm, prob)

                for i in range(len(self.states)):
                        for j in range(len(self.states)):
                            matrix[i][j] = elogproduct(matrix[i][j], -norm)

                xi_current.append(matrix)

            xi.append(xi_current)

        return xi

    def state_probs(self, alpha, beta):

        # probability of being at state i at time t
        gamma = []
        for seq in range(len(self.observations)):
            gamma_current = []
            for o in range(len(self.observations[seq])):
                norm = LOGZERO
                row = []
                for current in range(len(self.states)):
                        prob = elogproduct(alpha[seq][o][current],
                                           beta[seq][o][current])
                        row.append(prob)
                        norm = elogsum(norm, prob)
                for current in range(len(self.states)):
                        row[current] = elogproduct(row[current], -norm)

                gamma_current.append(row)

            gamma.append(gamma_current)

        return gamma

    def forward(self):

        # probability that the sequence from 0 to t was observed and
        # Markov process was at state j at time t
        alpha = []
        for seq in range(len(self.observations)):

            alpha_current = [[]]
            for s in range(len(self.states)):
                emissions = elog(self.emissions[s, self.observations[seq, 0]])
                init = elog(self.initial_p[s])
                alpha_current[0].append(elogproduct(emissions, init))
            for o in range(1, len(self.observations[seq])):
                row = []
                for current in range(len(self.states)):
                        total = LOGZERO
                        for prev in range(len(self.states)):
                            total = elogsum(
                                total,
                                elogproduct(alpha_current[o-1][prev],
                                            elog(self.transitions[prev][current])))
                        total = elogproduct(
                            total,
                            elog(self.emissions[current][self.observations[seq][o]]))
                        row.append(total)
                alpha_current.append(row)

            alpha.append(alpha_current)
        return alpha

    def backward(self):

        # probability that the sequence from t+1 to end was observed
        # and Markov process was at state j at time t
        beta = []
        for seq in range(len(self.observations)):
            beta_current = [[0] * len(self.states)
                            for i in range(len(self.observations[seq]))]
            for o in range(len(self.observations[seq]) - 2, -1, -1):
                for current in range(len(self.states)):
                        total = LOGZERO
                        for next in range(len(self.states)):
                            prod = elogproduct(
                                elog(self.emissions[next][self.observations[seq][o+1]]),
                                beta_current[o+1][next])
                            prod = elogproduct(elog(self.transitions[current][next]),
                                               prod)
                            total = elogsum(total, prod)
                        beta_current[o][current] = total

            beta.append(beta_current)

        return beta

    def calc_probs(self):
        obs_len = self.observations.size
        # T_states[i][j] gives state at position i - 1 that produces
        # maximum probability for path with state j at position i
        T_states = [[] for x in range(obs_len)]
        T_probs = [[] for x in range(obs_len)]

        # initialize
        for i in range(len(self.states)):
            T_probs[0].append(elogproduct(elog(self.initial_p[i]),
                                          elog(self.emissions[i, self.observations[0]])))

        # main loop of viterbi algorithm
        for pos in range(1, obs_len):
            for end_state in range(len(self.states)):
                max_prob = 'LOGZERO'
                max_state = -1
                for prev_state in range(len(self.states)):
                        trans_prob = self.transitions[prev_state][end_state]
                        emis_prob = self.emissions[end_state][self.observations[pos]]
                        prob = elogproduct(T_probs[pos - 1][prev_state],
                                           elogproduct(elog(trans_prob),
                                                       elog(emis_prob)))
                        if ecompare(prob, max_prob):
                            max_prob = prob
                            max_state = prev_state

                T_probs[pos].append(max_prob)
                T_states[pos].append(max_state)

        return T_probs, T_states

    def max_path(self, T_probs, T_states):

        max_path = []
        max_prob = 'LOGZERO'
        max_state = -1
        last_ind = len(self.observations) - 1
        # start by finding the state we are most probably in at the
        # last position
        for end_state in range(len(self.states)):
            p = T_probs[last_ind][end_state]
            if ecompare(p, max_prob):
                max_prob = p
                max_state = end_state
        max_path.append(max_state)

        # then retrace the most probable sequence of states that
        # preceded this state
        to_state = max_state
        # goes from last_ind to 1 backwards (inclusively)
        for i in range(last_ind, 0, -1):
            from_state = T_states[i][to_state]
            max_path.append(from_state)
            to_state = from_state

        max_path.reverse()
        return max_path

    def viterbi(self):
        self.observations = self.observations[0]
        T_probs, T_states = self.calc_probs()
        max_path = self.max_path(T_probs, T_states)
        return max_path

    def posterior_decoding(self):
        # probability of being in each state at each position, in
        # format of a list of lists
        alpha = self.forward()
        beta = self.backward()
        gamma = self.state_probs(alpha, beta)
        p = []
        for ind in range(len(self.observations)):
            p_ind = []
            for site in range(len(self.observations[ind])):
                p_ind_site = defaultdict(float)
                for state_ind in range(len(self.states)):
                        p_ind_site[self.states[state_ind]] = \
                            eexp(gamma[ind][site][state_ind])
                p_ind.append(p_ind_site)
            p.append(p_ind)
        return p
