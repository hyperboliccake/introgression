import hmm.hmm_bw as hmm
import pytest


@pytest.fixture
def hm3():
    hm = hmm.HMM()
    hm.set_hidden_states(list('ab'))
    hm.set_observed_states(list('NEA'))
    hm.set_emissions([{'N': 0.3, 'E': 0.6, 'A': 0.1},
                      {'N': 0.7, 'A': 0.2, 'E': 0.1}])
    hm.set_observations([list('NNENNENNEN'),
                         list('NNNANEENNN'),
                         list('NNNNNEENAN'),
                         list('NNENNEENEN')])
    hm.set_transitions([[0.5, 0.5], [0.3, 0.7]])
    hm.set_initial_p([0.2, 0.8])

    return hm


@pytest.fixture
def hm():
    hm = hmm.HMM()
    hm.set_hidden_states(['N', 'E'])
    hm.set_emissions([{'N': 0.3, 'E': 0.7},
                      {'N': 0.8, 'E': 0.2}])
    hm.set_observations([list('NNENNENNEN'),
                         list('NNNNNEENNN'),
                         list('NNENNEENEN')])
    hm.set_transitions([[0.5, 0.5], [0.3, 0.7]])
    hm.set_initial_p([0.2, 0.8])

    return hm
