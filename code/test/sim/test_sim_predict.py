from sim import sim_predict


def test_convert_predictions(hm):
    hm.observations = hm.observations
    vit = hm.viterbi()
    states = ['N', 'E']
    predicted = sim_predict.convert_predictions(vit, states)
    print(predicted)
    assert predicted == [states[i] for i in vit]
