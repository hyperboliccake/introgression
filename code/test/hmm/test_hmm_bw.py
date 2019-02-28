from hmm import hmm_bw as hmm
import pytest
from pytest import approx
import numpy as np


def test_init():
    hm = hmm.HMM()
    assert hm.hidden_states == []
    assert hm.transitions == []
    assert hm.emissions == []
    assert hm.observations == []
    assert hm.initial_p == []


def test_setters():
    hm = hmm.HMM()

    hm.set_hidden_states(['test'])
    assert hm.hidden_states == ['test']

    hm.set_transitions([[0, 1], [1, 0]])
    assert np.array_equal(hm.transitions, np.array([[0, 1], [1, 0]]))

    with pytest.raises(ValueError) as e:
        hm.set_transitions([[0, 0]])
    assert '[0, 0] 0' in str(e)

    hm.set_hidden_states(['1', '2'])
    hm.set_emissions([{'1': 1, '2': 0}, {'1': 0.5, '2': 0.5}])
    assert np.array_equal(hm.emissions, np.array([[1, 0], [0.5, 0.5]]))

    with pytest.raises(ValueError) as e:
        hm.set_emissions([{'1': 1, '2': 0}, {'1': 0, '2': 0}])
    assert '2 0' in str(e)

    hm.set_initial_p([0, 1, 0])
    assert np.array_equal(hm.initial_p, np.array([0, 1, 0]))

    with pytest.raises(AssertionError) as e:
        hm.set_initial_p([0, 0])
    assert '[0, 0] 0' in str(e)


def test_print_results(capsys, hm):
    hm.print_results(0, 1)
    captured = capsys.readouterr()
    out = captured.out.split('\n')

    assert out[0] == 'Iterations: 0'
    assert out[2] == 'Log Likelihood:'
    assert out[5] == 'Initial State Probabilities:'
    assert out[9] == 'Transition Probabilities:'
    assert out[15] == 'Emission Probabilities:'

    assert float(out[3]) == 1

    assert float(out[6].split('=')[1]) == 0.2
    assert float(out[7].split('=')[1]) == 0.8

    assert float(out[10].split('=')[1]) == 0.5
    assert float(out[11].split('=')[1]) == 0.5
    assert float(out[12].split('=')[1]) == 0.3
    assert float(out[13].split('=')[1]) == 0.7

    assert float(out[16].split('=')[1]) == 0.7
    assert float(out[17].split('=')[1]) == 0.3
    assert float(out[18].split('=')[1]) == 0.2
    assert float(out[19].split('=')[1]) == 0.8


def test_go(capsys, hm):
    hm.go()
    # get output from last report
    out = capsys.readouterr().out.split('\n')[-23:]
    assert out[0] == 'Iterations: 2'
    assert out[2] == 'Log Likelihood:'
    assert out[5] == 'Initial State Probabilities:'
    assert out[9] == 'Transition Probabilities:'
    assert out[15] == 'Emission Probabilities:'

    assert round(float(out[3]), 1) == -17.8

    assert round(float(out[6].split('=')[1]), 3) == 0.033
    assert round(float(out[7].split('=')[1]), 2) == 0.97

    assert round(float(out[10].split('=')[1]), 2) == 0.41
    assert round(float(out[11].split('=')[1]), 2) == 0.59
    assert round(float(out[12].split('=')[1]), 2) == 0.32
    assert round(float(out[13].split('=')[1]), 2) == 0.68

    assert round(float(out[16].split('=')[1]), 2) == 0.62
    assert round(float(out[17].split('=')[1]), 2) == 0.38
    assert round(float(out[18].split('=')[1]), 2) == 0.15
    assert round(float(out[19].split('=')[1]), 2) == 0.85


def test_intial_likelihood(hm):
    alpha = hm.forward()
    LL = hm.log_likelihood(alpha)

    iter_LL = 0
    for seq in range(len(hm.observations)):
        LL_seq = np.copy(alpha[seq, -1, 0])
        for i in range(1, len(hm.hidden_states)):
            LL_seq = np.logaddexp(LL_seq, alpha[seq, -1, i])
        iter_LL += LL_seq

    assert LL == approx(iter_LL)


def test_initial_probabilities(hm):
    gamma = hm.state_probs(hm.forward(), hm.backward())
    pi = hm.initial_probabilities(gamma)
    iter_pi = np.copy(gamma[0, 0, :])

    for i in range(len(hm.hidden_states)):
        for seq in range(1, len(hm.observations)):
            iter_pi[i] = np.logaddexp(iter_pi[i], gamma[seq, 0, i])
    iter_pi = [np.exp(x) / len(hm.observations) for x in iter_pi]

    assert pi == approx(iter_pi)


def test_transition_probabilities(hm):
    alpha = hm.forward()
    beta = hm.backward()
    gamma = hm.state_probs(alpha, beta)
    xi = hm.bw(alpha, beta)

    trans = hm.transition_probabilities(xi, gamma)

    iter_trans = []
    for i in range(len(hm.hidden_states)):
        row = []
        for j in range(len(hm.hidden_states)):
            num = np.NINF
            den = np.NINF
            for seq in range(len(hm.observations)):
                num_seq = xi[seq][0][i][j]
                den_seq = gamma[seq][0][i]
                for o in range(1, len(hm.observations[seq]) - 1):
                    # xi is probability of probability
                    # of being at state i at time t and
                    # state j at time t+1
                    num_seq = np.logaddexp(num_seq, xi[seq][o][i][j])
                    # gamma is probability of being in
                    # state i at time t
                    den_seq = np.logaddexp(den_seq, gamma[seq][o][i])
                # add the current sequence contribution to total
                num = np.logaddexp(num, num_seq)
                den = np.logaddexp(den, den_seq)
            row.append(np.exp(num - den))
        iter_trans.append(row)

    assert iter_trans == approx(trans)


def test_emission_probabilities(hm3):
    hm = hm3
    gamma = hm.state_probs(hm.forward(), hm.backward())
    em = hm.emission_probabilities(gamma)

    # emission probabilities
    iter_em = []
    for state in range(len(hm.hidden_states)):
        d = []
        for symbol in range(len(hm.emissions[state])):
            num = np.NINF
            den = np.NINF
            for seq in range(len(hm.observations)):
                num_seq = np.NINF
                den_seq = np.NINF
                for o in range(len(hm.observations[seq])):
                    if hm.observations[seq][o] == symbol:
                        # gamma is probability of being in state i at time t
                        num_seq = np.logaddexp(num_seq, gamma[seq][o][state])
                    den_seq = np.logaddexp(den_seq, gamma[seq][o][state])

                # add the current sequence contribution to total
                num = np.logaddexp(num, num_seq)
                den = np.logaddexp(den, den_seq)
            d.append(np.exp(num - den))
        iter_em.append(d)

    iter_em = np.array(iter_em)
    assert em == approx(iter_em)

    den = np.logaddexp.reduce(gamma, axis=0)
    den = np.logaddexp.reduce(den, axis=0)

    obs = np.array([i == hm.observations for i in range(len(hm.observed_states))])
    obs = np.moveaxis(obs, [0, 1, 2], [2, 0, 1])
    gam = np.where(obs[:, :, None, :], gamma[:, :, :, None], np.NINF)
    num = np.logaddexp.reduce(np.logaddexp.reduce(gam))

    assert em == approx(np.exp(num - den[:, None]))


def test_forward(hm):
    alpha = hm.forward()

    alpha_iter = []
    for seq in range(len(hm.observations)):

        alpha_current = [[]]
        for s in range(len(hm.hidden_states)):
            alpha_current[0].append(np.log(
                hm.emissions[s][hm.observations[seq][0]]
                * hm.initial_p[s]))
        for o in range(1, len(hm.observations[seq])):
            row = []
            for current in range(len(hm.hidden_states)):
                    total = -np.inf
                    for prev in range(len(hm.hidden_states)):
                        total = np.logaddexp(
                            total,
                            alpha_current[o-1][prev] +
                            np.log(hm.transitions[prev][current]))
                    total += np.log(
                        hm.emissions[current][hm.observations[seq][o]])
                    row.append(total)
            alpha_current.append(row)

        alpha_iter.append(alpha_current)
    alpha_iter = np.array(alpha_iter)

    assert alpha == approx(alpha_iter)


def test_backward(hm):
    beta = hm.backward()

    beta_iter = []
    for seq in range(len(hm.observations)):
        beta_current = [[0] * len(hm.hidden_states)
                        for i in range(len(hm.observations[seq]))]
        for o in range(len(hm.observations[seq]) - 2, -1, -1):
            for current in range(len(hm.hidden_states)):
                total = np.NINF
                for next in range(len(hm.hidden_states)):
                    prod = np.log(
                        hm.emissions[next][hm.observations[seq][o+1]]) +\
                        beta_current[o+1][next]
                    prod += np.log(hm.transitions[current][next])
                    total = np.logaddexp(total, prod)
                beta_current[o][current] = total

        beta_iter.append(beta_current)
    beta_iter = np.array(beta_iter)

    assert beta == approx(beta_iter)


def test_state_probs(hm):
    alpha = hm.forward()
    beta = hm.backward()
    gamma = hm.state_probs(alpha, beta)

    iter_gamma = []
    for seq in range(len(hm.observations)):
        gamma_current = []
        for o in range(len(hm.observations[seq])):
            norm = np.NINF
            row = []
            for current in range(len(hm.hidden_states)):
                prob = alpha[seq][o][current] + beta[seq][o][current]
                row.append(prob)
                norm = np.logaddexp(norm, prob)
            for current in range(len(hm.hidden_states)):
                row[current] = row[current] - norm

            gamma_current.append(row)

        iter_gamma.append(gamma_current)

    iter_gamma = np.array(iter_gamma)

    assert gamma == approx(iter_gamma)


def test_bw(hm):
    alpha = hm.forward()
    beta = hm.backward()
    xi = hm.bw(alpha, beta)

    xi_iter = []
    for seq in range(len(hm.observations)):
        xi_current = []
        for o in range(len(hm.observations[seq]) - 1):
            norm = np.NINF
            matrix = [[np.NINF] * len(hm.hidden_states)
                      for i in range(len(hm.hidden_states))]
            for i in range(len(hm.hidden_states)):
                for j in range(len(hm.hidden_states)):
                    prob = np.log(
                        hm.emissions[j][hm.observations[seq][o+1]]
                    ) + beta[seq][o+1][j]
                    prob = np.log(hm.transitions[i][j]) + prob
                    prob = alpha[seq][o][i] + prob
                    matrix[i][j] = prob
                    norm = np.logaddexp(norm, prob)

            for i in range(len(hm.hidden_states)):
                    for j in range(len(hm.hidden_states)):
                        matrix[i][j] = matrix[i][j] - norm

            xi_current.append(matrix)

        xi_iter.append(xi_current)

    xi_iter = np.array(xi_iter)
    xi_iter[np.isnan(xi_iter)] = np.NINF

    assert xi == approx(xi_iter)


# TODO make function call take a sequence to consider
def test_calc_probs(hm):
    hm.observations = hm.observations[0]
    probs, states = hm.calculate_max_states()

    obs_len = hm.observations.size
    iter_states = [[] for x in range(obs_len)]
    iter_probs = [[] for x in range(obs_len)]

    # initialize
    for i in range(len(hm.hidden_states)):
        iter_probs[0].append(np.log(hm.initial_p[i]) +
                             np.log(hm.emissions[i, hm.observations[0]]))
        iter_states[0].append(-1)

    # main loop of viterbi algorithm
    for pos in range(1, obs_len):
        for end_state in range(len(hm.hidden_states)):
            max_prob = np.NINF
            max_state = -1
            for prev_state in range(len(hm.hidden_states)):
                    trans_prob = hm.transitions[prev_state][end_state]
                    emis_prob = hm.emissions[end_state][hm.observations[pos]]
                    prob = iter_probs[pos - 1][prev_state] + \
                        np.log(trans_prob) + np.log(emis_prob)
                    if prob > max_prob:
                        max_prob = prob
                        max_state = prev_state

            iter_probs[pos].append(max_prob)
            iter_states[pos].append(max_state)

    iter_probs = np.array(iter_probs)
    iter_states = np.array(iter_states)

    assert probs == approx(iter_probs)
    assert states == approx(iter_states)


def test_max_path(hm):
    hm.observations = hm.observations[0]
    prob, states = hm.calculate_max_states()
    path = hm.max_path(prob, states)

    max_path = []
    max_prob = np.NINF
    max_state = -1
    last_ind = len(hm.observations) - 1
    # start by finding the state we are most probably in at the
    # last position
    for end_state in range(len(hm.hidden_states)):
        p = prob[last_ind][end_state]
        if p > max_prob:
            max_prob = p
            max_state = end_state
    max_path.append(max_state)

    # then retrace the most probable sequence of states that
    # preceded this state
    to_state = max_state
    # goes from last_ind to 1 backwards (inclusively)
    for i in range(last_ind, 0, -1):
        from_state = states[i][to_state]
        max_path.append(from_state)
        to_state = from_state

    max_path.reverse()
    assert path == approx(max_path)

    hm.observations = [hm.observations]
    assert np.array(hm.viterbi()) == approx(max_path)


def test_posterior_decoding(hm):
    post = hm.posterior_decoding()

    gamma = np.exp(hm.state_probs(hm.forward(), hm.backward()))

    assert post == approx(gamma)
