from sim import sim_predict
from hmm import hmm_bw as hmm
import pytest


@pytest.fixture
def hm():
    hm = hmm.HMM()
    hm.set_states(['N', 'E'])
    hm.set_emissions([{'N': 0.3, 'E': 0.7},
                      {'N': 0.8, 'E': 0.2}])
    hm.set_observations([list('NNENNENNEN'),
                         list('NNNNNEENNN'),
                         list('NNENNEENEN')])
    hm.set_transitions([[0.5, 0.5], [0.3, 0.7]])
    hm.set_initial_p([0.2, 0.8])

    return hm


def test_convert_predictions(hm):
    hm.observations = hm.observations
    vit = hm.viterbi()
    states = ['N', 'E']
    predicted = sim_predict.convert_predictions(vit, states)
    print(predicted)
    assert predicted == [states[i] for i in vit]
