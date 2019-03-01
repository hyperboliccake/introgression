import hmm.hmm_bw as hmm
import math
import random
import pytest
from pytest import approx
import numpy as np


def test_elog():
    assert hmm.elog(0) == hmm.LOGZERO
    assert hmm.elog(hmm.LOGZERO) == hmm.LOGZERO
    random.seed(123)
    for i in range(10):
        num = random.random()
        assert hmm.elog(num) == math.log(num)


def test_elogproduct():
    assert hmm.elogproduct(0, hmm.LOGZERO) == hmm.LOGZERO
    assert hmm.elogproduct(hmm.LOGZERO, 0) == hmm.LOGZERO
    random.seed(123)
    for i in range(10):
        m = random.random()
        n = random.random()
        assert hmm.elogproduct(m, n) == m + n


def test_elogsum():
    assert hmm.elogsum(0, hmm.LOGZERO) == 0
    assert hmm.elogsum(hmm.LOGZERO, 0) == 0

    assert hmm.elogsum(1, 1) == 1 + math.log(2)

    random.seed(123)
    for i in range(10):
        m = random.random()
        n = random.random()
        assert hmm.elogsum(m, n) ==\
            math.log(1 + math.exp(-1*abs(m-n))) + max(m, n)


def test_eexp():
    assert hmm.eexp(hmm.LOGZERO) == 0

    random.seed(123)
    for i in range(10):
        m = random.random()
        assert hmm.eexp(m) == math.exp(m)


def test_ecompare():
    assert hmm.ecompare(0, hmm.LOGZERO) is True
    assert hmm.ecompare(hmm.LOGZERO, hmm.LOGZERO) is False
    assert hmm.ecompare(hmm.LOGZERO, 0) is False

    for i in range(10):
        m = random.random()
        n = random.random()
        assert hmm.ecompare(m, n) == (m > n)


def test_isclose():
    assert True  # TODO replace with math.isclose in python 3


def test_init():
    hm = hmm.HMM()
    assert hm.states == []
    assert hm.trans == []
    assert hm.emis == []
    assert hm.obs == []
    assert hm.init == []


def test_setters():
    hm = hmm.HMM()

    # TODO probably need better type checking here
    hm.set_states('test')
    assert hm.states == 'test'

    hm.set_obs('test')
    assert hm.obs == 'test'

    hm.set_trans([[0, 1], [1, 0]])
    assert hm.trans == [[0, 1], [1, 0]]

    with pytest.raises(AssertionError) as e:
        hm.set_trans([[0, 0]])
    assert '[0, 0] 0' in str(e)

    hm.set_emis([{'1': 1, '2': 0}, {'1': 0, '2': 1}])
    assert hm.emis == [{'1': 1, '2': 0}, {'1': 0, '2': 1}]

    with pytest.raises(AssertionError) as e:
        hm.set_emis([[0, 0]])
    assert '[0, 0] 0' in str(e)

    hm.set_init([0, 1, 0])
    assert hm.init == [0, 1, 0]

    with pytest.raises(AssertionError) as e:
        hm.set_init([0, 0])
    assert '[0, 0] 0' in str(e)


@pytest.fixture
def hm():
    hm = hmm.HMM()
    hm.set_states(['N', 'E'])
    hm.set_obs([list('NNENNENNEN'),
                list('NNNNNEENNN'),
                list('NNENNEENEN')])
    hm.set_trans([[0.5, 0.5], [0.3, 0.7]])
    hm.set_emis([{'N': 0.3, 'E': 0.7},
                 {'N': 0.8, 'E': 0.2}])
    hm.set_init([0.2, 0.8])

    return hm


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


def test_forward(hm):
    for seqnum in range(len(hm.obs)):
        alpha = hm.forward()[seqnum]
        seq = hm.obs[seqnum]

        # emis * init for state 1, seq N
        assert alpha[0][0] == approx(math.log(0.3 * 0.2))
        assert alpha[0][1] == approx(math.log(0.8 * 0.8))

        for i in range(1, len(alpha)):
            alpha_a = np.dot(np.exp(alpha[i-1]), hm.trans)
            assert alpha[i][0] == approx(
                math.log(hm.emis[0][seq[i]] * alpha_a[0]))
            assert alpha[i][1] == approx(
                math.log(hm.emis[1][seq[i]] * alpha_a[1]))


def test_backward(hm):
    for seqnum in range(len(hm.obs)):
        beta = hm.backward()[seqnum]
        # reverse so recursion is easier
        beta.reverse()
        seq = hm.obs[seqnum]
        seq.reverse()

        assert beta[0][0] == math.log(1)
        assert beta[0][0] == math.log(1)

        for i in range(1, len(beta)):
            beta_ab = np.dot(hm.trans,
                             np.multiply(np.exp(beta[i-1]),
                                         [hm.emis[j][seq[i-1]]
                                          for j in range(2)]))
            assert beta[i][0] == approx(math.log(beta_ab[0]))


def test_state_probs(hm):
    alpha = hm.forward()
    beta = hm.backward()
    gamma = hm.state_probs(alpha, beta)

    alpha = np.exp(alpha)
    beta = np.exp(beta)
    ab = np.multiply(alpha, beta)

    assert gamma == approx(np.log(
        ab / np.sum(ab, axis=2)[:, :, np.newaxis]))


def test_bw(hm):
    Alpha = hm.forward()
    Beta = hm.backward()
    Xi = hm.bw(Alpha, Beta)

    for seqnum in range(len(hm.obs)):
        # offset to match time t with t+1 in beta and seq
        alpha = np.exp(Alpha[seqnum][:-1])
        beta = np.exp(Beta[seqnum][1:])
        seq = hm.obs[seqnum][1:]
        xi = np.exp(Xi[seqnum])
        # make emission a dir of lists instead of list of dirs
        emis = {i: [hm.emis[j][i] for j in range(2)]
                for i in hm.states}

        # build array of emission based on sequence
        emis = np.array([emis[n] for n in seq])

        # element-wise multiplication with weird broadcasting
        num = alpha[:, :, None] * beta[:, None, :]\
            * hm.trans * emis[:, None, :]

        # normalize by summing along axis 1 and 2, broadcast division
        assert xi == approx(num / np.sum(num, axis=(1, 2))[:, None, None])


# TODO make function call take a sequence to consider
def test_calc_probs(hm):
    # seems like this is only supported for a single sequence
    hm.obs = hm.obs[0]
    prob, states = hm.calc_probs()
    prob = np.exp(prob)

    # make emission a dir of lists instead of list of dirs
    emis = {i: [hm.emis[j][i] for j in range(2)]
            for i in hm.states}

    # build array of emission based on sequence
    emis = np.array([emis[n] for n in hm.obs])

    p = [list(hm.init * emis[0])]
    trans_emis = np.array(hm.trans)[None, :, :] * emis[:, None, :]
    s = [[]]

    for i in range(1, len(emis)):
        probs = p[-1] * np.transpose(trans_emis[i])
        p.append(list(np.max(probs, axis=1)))
        s.append(list(np.argmax(probs, axis=1)))

    assert states == s
    for i in range(len(prob)):
        assert prob[i] == approx(p[i])


def test_max_path(hm):
    hm.obs = hm.obs[0]
    prob, states = hm.calc_probs()
    path = hm.max_path(prob, states)

    p = [np.argmax(prob.pop())]
    states = states[1:]
    while states:
        p.append(states.pop()[p[-1]])

    p.reverse()

    assert p == path

    assert hm.viterbi() == p


def test_posterior_decoding(hm):
    post = hm.posterior_decoding()

    gamma = np.exp(hm.state_probs(hm.forward(), hm.backward()))

    from collections import defaultdict
    pdd = [[defaultdict(float,
                        {key: value for key, value in zip(hm.states, entry)})
            for entry in seq] for seq in gamma]

    for i in range(len(post)):
        for j in range(len(post[i])):
            for k in post[i][j].keys():
                assert post[i][j][k] == approx(pdd[i][j][k])
