    


def predict_introgressed_hmm_old(seqs_filled, third):

    # create a hidden markov model and determine which reference genome we
    # are most likely to be in at each variant site
    hmm = HMM()
        
    hmm.set_obs(seqs_filled)

    # only cer and par
    if third == 'none':
        hmm.set_init([.5, .5])
        hmm.set_states(['cer', 'par'])

        hmm_trans_cer_par = .0005

        hmm.set_trans([\
                [1-hmm_trans_cer_par,hmm_trans_cer_par],\
                [hmm_trans_cer_par,1-hmm_trans_cer_par],\
                    ])

        emis_cer = {\
            hmm_symbol['cer']:.5,\
                hmm_symbol['cerpar']:.4998,\
                hmm_symbol['par']:.0001,\
                hmm_symbol['none']:.0001\
                }
        assert sum(emis_cer.values()) == 1, sum(emis_cer.values())

        emis_par = {\
            hmm_symbol['par']:.5,\
                hmm_symbol['cerpar']:.4998,\
                hmm_symbol['cer']:.0001,\
                hmm_symbol['none']:.0001\
                }
        assert sum(emis_par.values()) == 1, sum(emis_par.values())
    
        hmm.set_emis([emis_cer, emis_par])
    
    elif third == 'bay':

        hmm.set_init([1/3.,1/3.,1/3.])
        # anything (of the correct length) would work here, because really
        # we just use indices; difference between bayanus and unknown is
        # in whether it matches bayanus sequence or just doesn't match
        # cerevisiae and paradoxus
        hmm.set_states(['cer', 'par', 'bay'])
        # a little weird because we know recombination rate we used
        hmm_trans_cer_par = .0005
        hmm_trans_cer_bay = .0005
        hmm_trans_par_bay = .0005
    
        # hmm.set_trans({'cer':{'cer':1-hmm_trans, 'par':hmm_trans},
        # 'par':{'cer':hmm_trans, 'par':1-hmm_trans}})
        hmm.set_trans([\
                [1-(hmm_trans_cer_par-hmm_trans_cer_bay),hmm_trans_cer_par,hmm_trans_cer_bay],\
                    [hmm_trans_cer_par,1-(hmm_trans_cer_par-hmm_trans_par_bay),hmm_trans_par_bay],\
                    [hmm_trans_cer_bay,hmm_trans_par_bay,1-(hmm_trans_cer_bay-hmm_trans_par_bay)]\
                    ])
        
        # 0    1       2       3          4    5       6    7
        # cer, cerpar, cerbay, cerparbay, par, parbay, bay, none
        emis_cer = {\
            hmm_symbol['cer']:.5,\
                hmm_symbol['cerpar']:.25,\
                hmm_symbol['cerbay']:.15,\
                hmm_symbol['cerparbay']:.0995,\
                hmm_symbol['par']:.0001,\
                hmm_symbol['parbay']:.0001,\
                hmm_symbol['bay']:.0001,\
                hmm_symbol['none']:.0001\
                }
        assert sum(emis_cer.values()) == 1, sum(emis_cer.values())
        
        emis_par = {\
            hmm_symbol['par']:.5,\
                hmm_symbol['cerpar']:.25,\
                hmm_symbol['parbay']:.15,\
                hmm_symbol['cerparbay']:.0995,\
                hmm_symbol['cer']:.0001,\
                hmm_symbol['cerbay']:.0001,\
                hmm_symbol['bay']:.0001,\
                hmm_symbol['none']:.0001\
                }
        assert sum(emis_par.values()) == 1, sum(emis_par.values())

        emis_bay = {\
            hmm_symbol['bay']:.5,\
                hmm_symbol['parbay']:.25,\
                hmm_symbol['cerbay']:.15,\
                hmm_symbol['cerparbay']:.0995,\
                hmm_symbol['par']:.0001,\
                hmm_symbol['cerpar']:.0001,\
                hmm_symbol['cer']:.0001,\
                hmm_symbol['none']:.0001\
                }
        assert sum(emis_bay.values()) == 1, sum(emis_bay.values())

        hmm.set_emis([emis_cer, emis_par, emis_bay])

    elif third == 'unknown':

        hmm.set_init([1/3.,1/3.,1/3.])
        hmm.set_states(['cer', 'par', 'unk'])

        hmm_trans_cer_par = .0005
        hmm_trans_cer_unk = .0005
        hmm_trans_par_unk = .0005

        hmm.set_trans([\
                [1-(hmm_trans_cer_par-hmm_trans_cer_unk),hmm_trans_cer_par,hmm_trans_cer_unk],\
                    [hmm_trans_cer_par,1-(hmm_trans_cer_par-hmm_trans_par_unk),hmm_trans_par_unk],\
                    [hmm_trans_cer_unk,hmm_trans_par_unk,1-(hmm_trans_cer_unk-hmm_trans_par_unk)]\
                    ])

        emis_cer = {\
            hmm_symbol['cer']:.5,\
                hmm_symbol['cerpar']:.4998,\
                hmm_symbol['par']:.0001,\
                hmm_symbol['none']:.0001\
                }
        assert sum(emis_cer.values()) == 1, sum(emis_cer.values())
        
        emis_par = {\
            hmm_symbol['par']:.5,\
                hmm_symbol['cerpar']:.4998,\
                hmm_symbol['cer']:.0001,\
                hmm_symbol['none']:.0001\
                }
        assert sum(emis_par.values()) == 1, sum(emis_par.values())

        emis_unk = {\
            hmm_symbol['par']:.0001,\
                hmm_symbol['cerpar']:.0001,\
                hmm_symbol['cer']:.0001,\
                hmm_symbol['none']:.9997\
                }
        assert sum(emis_unk.values()) == 1, sum(emis_unk.values())

        hmm.set_emis([emis_cer, emis_par, emis_unk])

    else:
        print 'option for third incorrect'
        sys.exit()
                
    # Baum-Welch parameter estimation
    hmm.go()

    predicted = []
    for i in range(len(seqs_filled)):
        obs_seq = seqs_filled[i]
        hmm.set_obs(obs_seq)
        predicted.append(hmm.viterbi())

    return predicted, hmm
