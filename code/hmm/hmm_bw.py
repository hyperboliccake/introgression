import math
import sys

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

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

class HMM:

     def __init__(self):

          self.states = []
          self.trans = []
          self.emis = []
          self.obs = []
          self.init = []

     def set_states(self, states):

          self.states = states

     def set_trans(self, trans):
         
          self.trans = trans
          for x in self.trans:
               t = 0
               for i in x:
                    t += i
               assert isclose(t, 1), str(x) + ' ' + str(t)

     def set_emis(self, emis):
         
          self.emis = emis
          for x in self.emis:
               t = 0
               for i in x:
                    t += x[i]
               assert isclose(t, 1), str(x) + ' ' + str(t)

     def set_obs(self, obs):
         
          # one sequence for Viterbi, or list of observed sequences
          # for Baum-Welch; TODO maybe change this weirdness?
          self.obs = obs

     def set_init(self, init):
          
          self.init = init
          assert isclose(sum(self.init), 1), str(self.init) + ' ' + str(sum(self.init))

     def print_results(self, num_its, LL):

          print "Iterations:", num_its
          print
          print "Log Likelihood:"
          print '%.30e' % LL
          print
          print "Initial State Probabilities:"
          for i in xrange(len(self.states)):
               print str(self.states[i]) + '=' + '%.30e' % self.init[i]
          print
          print "Transition Probabilities:"
          for i in xrange(len(self.states)):
               for j in xrange(len(self.states)):
                    print str(self.states[i]) + ',' + str(self.states[j]) + '=' + '%.30e' % self.trans[i][j]
          print          
          print "Emission Probabilities:"
          for i in xrange(len(self.states)):
               keys = sorted(self.emis[i].keys())
               for k in keys:
                    print str(self.states[i]) + ',' + k + '=' + '%.30e' % self.emis[i][k]
          print

          sys.stdout.flush()

     def go(self, improvement_frac = .01, num_its = 0, prev_LL = 0):

          # calculate current log likelihood
          print "calculating alpha"
          alpha = self.forward()

          # for multiple observations, product of LL
          # (note not base 2 log)
          LL_all = []
          LL = math.log(1) # because we're taking product
          for seq in range(len(self.obs)):
               LL_seq = LOGZERO # because we're taking sum
               for i in xrange(len(self.states)):
                    LL_seq = elogsum(LL_seq, alpha[seq][-1][i])
               LL_all.append(LL_seq)
               LL = elogproduct(LL, LL_seq) 
          self.print_results(num_its, LL)

          # continue until log likelihood has stopped increasing much
          threshold = improvement_frac * abs(LL)
          while num_its < 1 or LL - prev_LL > threshold:

               print 'Iteration',  num_its
                    
               print "calculating beta"
               beta = self.backward()
               print "calculating gamma"
               gamma = self.state_probs(alpha, beta)
               print "calculating xi"
               xi = self.bw(alpha, beta)

               print "updating parameters"

               # initial probabilities
               pi = [LOGZERO for i in xrange(len(self.states))]
               for i in xrange(len(self.states)):
                    for seq in xrange(len(self.obs)):
                         pi[i] = elogsum(pi[i], gamma[seq][0][i])
               pi = [eexp(x) / len(self.obs) for x in pi]
               """
               pi = [0. for i in xrange(len(self.states))]
               for seq in xrange(len(self.obs)):
                    pi_seq = []
                    for i in xrange(len(self.states)):
                         pi_seq.append(gamma[seq][0][i])
                    pi_seq = [elogproduct(x, -LL_all[seq]) for x in pi_seq]
                    pi = [elogsum(pi[i], pi_seq[i]) for i in range(len(self.states))]
               sum_pi = reduce(elogsum, pi)
               pi = [eexp(elogproduct(x, -sum_pi)) for x in pi]
               """
               # transition probabilities
               a = []
               for i in xrange(len(self.states)):
                    row = []
                    for j in xrange(len(self.states)):
                         num = LOGZERO
                         den = LOGZERO
                         for seq in range(len(self.obs)):
                              num_seq = LOGZERO
                              den_seq = LOGZERO
                              for o in xrange(len(self.obs[seq]) - 1):
                                   # xi is probability of probability
                                   # of being at state i at time t and
                                   # state j at time t+1
                                   num_seq = elogsum(num_seq, xi[seq][o][i][j])
                                   # gamma is probability of being in
                                   # state i at time t
                                   den_seq = elogsum(den_seq, gamma[seq][o][i])
                              # weight numerator and denominator for the
                              # current observation sequence by 1 / P(seq | model)
                              #num_seq = elogproduct(num_seq, -LL_all[seq])
                              #den_seq = elogproduct(den_seq, -LL_all[seq])
                              # add the current sequence contribution to total
                              num = elogsum(num, num_seq)
                              den = elogsum(den, den_seq)
                         assert den != LOGZERO, \
                             'probably something wrong with initial parameter values'
                         row.append(eexp(elogproduct(num, -den)))
                    a.append(row)
     
               # emission probabilities
               b = []
               for state in xrange(len(self.states)):
                    d = {}
                    for symbol in self.emis[state].keys():
                         num = LOGZERO
                         den = LOGZERO
                         for seq in range(len(self.obs)):
                              num_seq = LOGZERO
                              den_seq = LOGZERO
                              for o in xrange(len(self.obs[seq])):
                                   if self.obs[seq][o] == symbol:
                                        # gamma is probability of
                                        # being in state i at time t
                                        num_seq = elogsum(num_seq, gamma[seq][o][state])
                                   den_seq = elogsum(den_seq, gamma[seq][o][state])
                              # weight numerator and denominator for the
                              # current observation sequence by 1 / P(seq | model)
                              #num_seq = elogproduct(num_seq, -LL_all[seq])
                              #den_seq = elogproduct(den_seq, -LL_all[seq])
                              # add the current sequence contribution to total
                              num = elogsum(num, num_seq)
                              den = elogsum(den, den_seq)
                         assert den != LOGZERO, \
                             'probably something wrong with initial parameter values'
                         d[symbol] = eexp(elogproduct(num, -den))
                    b.append(d)
          
               self.init = pi
               self.trans = a
               self.emis = b

               assert isclose(sum(self.init), 1), str(sum(self.init)) + ' ' + str(self.init)
               for t in self.trans:
                    assert isclose(sum(t), 1), str(sum(t)) + ' ' + str(t)
               for e in self.emis:
                    assert isclose(sum(e.values()), 1), str(sum(e.values())) + ' ' + str(e)


               num_its += 1

               print "calculating alpha"
               alpha = self.forward()
               
               prev_LL = LL

               LL_all = []
               LL = math.log(1) # because we're taking product
               for seq in range(len(self.obs)):
                    LL_seq = LOGZERO # because we're taking sum
                    for i in xrange(len(self.states)):
                         LL_seq = elogsum(LL_seq, alpha[seq][-1][i])
                    LL_all.append(LL_seq)
                    LL = elogproduct(LL, LL_seq) 

               # print results for every iteration
               self.print_results(num_its, LL)

               assert LL > prev_LL or isclose(LL, prev_LL), \
                    str(LL) + ' ' + str(prev_LL)

               
          print "finished in", num_its, 'iterations'

     def bw(self, alpha, beta):

          # probability of being at state i at time t and state j at
          # time t+1
          xi = []
          for seq in range(len(self.obs)):
               xi_current = []
               for o in xrange(len(self.obs[seq]) - 1):
                    norm = LOGZERO
                    matrix = [[LOGZERO] * len(self.states) for i in xrange(len(self.states))]
                    for i in xrange(len(self.states)):
                         for j in xrange(len(self.states)):
                              #emis_prob = -1
                              #if self.obs[seq][o+1] == self.state_seqs[j][o+1]
                              prob = elogproduct(elog(self.emis[j][self.obs[seq][o+1]]), beta[seq][o+1][j])
                              prob = elogproduct(elog(self.trans[i][j]), prob)
                              prob = elogproduct(alpha[seq][o][i], prob)
                              matrix[i][j] = prob
                              norm = elogsum(norm, prob)

                    for i in xrange(len(self.states)):
                         for j in xrange(len(self.states)):
                              matrix[i][j] = elogproduct(matrix[i][j], -norm)
                    
                    xi_current.append(matrix)

               xi.append(xi_current)

          return xi

     def state_probs(self, alpha, beta):

          # probability of being at state i at time t
          gamma = []
          for seq in range(len(self.obs)):
               gamma_current = []
               for o in xrange(len(self.obs[seq])):
                    norm = LOGZERO
                    row = []
                    for current in xrange(len(self.states)):
                         prob = elogproduct(alpha[seq][o][current], beta[seq][o][current])
                         row.append(prob)
                         norm = elogsum(norm, prob)
                    for current in xrange(len(self.states)):
                         row[current] = elogproduct(row[current], -norm)

                    gamma_current.append(row)

               gamma.append(gamma_current)

          return gamma

     def forward(self):

          # probability that the sequence from 0 to t was observed and
          # Markov process was at state j at time t
          alpha = []
          for seq in range(len(self.obs)):

               alpha_current = [[]]
               for s in xrange(len(self.states)):
                    emis = elog(self.emis[s][self.obs[seq][0]])
                    init = elog(self.init[s])
                    alpha_current[0].append(elogproduct(emis, init))
               for o in xrange(1, len(self.obs[seq])):
                    row = []
                    for current in xrange(len(self.states)):
                         total = LOGZERO
                         for prev in xrange(len(self.states)):
                              total = elogsum(total, elogproduct(alpha_current[o-1][prev], elog(self.trans[prev][current])))
                         total = elogproduct(total, elog(self.emis[current][self.obs[seq][o]]))
                         row.append(total)
                    alpha_current.append(row)

               alpha.append(alpha_current)
          return alpha

     def backward(self):

          # probability that the sequence from t+1 to end was observed
          # and Markov process was at state j at time t
          beta = []
          for seq in range(len(self.obs)):
               beta_current = [[0] * len(self.states) for i in xrange(len(self.obs[seq]))]
               for o in xrange(len(self.obs[seq]) - 2, -1, -1):
                    for current in xrange(len(self.states)):
                         total = LOGZERO
                         for next in xrange(len(self.states)):
                              prod = elogproduct(elog(self.emis[next][self.obs[seq][o+1]]), beta_current[o+1][next])
                              prod = elogproduct(elog(self.trans[current][next]), prod)
                              total = elogsum(total, prod)
                         beta_current[o][current] = total

               beta.append(beta_current)

          return beta

     def calc_probs(self):

          obs_len = len(self.obs)
          # T_states[i][j] gives state at position i - 1 that produces
          # maximum probability for path with state j at position i
          T_states = [[] for x in xrange(obs_len)]
          T_probs = [[] for x in xrange(obs_len)]

          # initialize
          for i in range(len(self.states)):
               T_probs[0].append(elogproduct(elog(self.init[i]), elog(self.emis[i][self.obs[0]])))
               
          # main loop of viterbi algorithm
          for pos in xrange(1, obs_len):
               #if pos % 10000 == 0:
               #     print pos
               for end_state in range(len(self.states)):
                    max_prob = 'LOGZERO'
                    max_state = -1
                    for prev_state in range(len(self.states)):
                         trans_prob = self.trans[prev_state][end_state]
                         emis_prob = self.emis[end_state][self.obs[pos]]
                         prob = elogproduct(T_probs[pos - 1][prev_state],\
                                                 elogproduct(elog(trans_prob), elog(emis_prob)))
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
          last_ind = len(self.obs) - 1
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
          for i in xrange(last_ind, 0, -1):
               from_state = T_states[i][to_state]
               max_path.append(from_state)
               to_state = from_state

          max_path.reverse()
          return max_path

     def viterbi(self):
          
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
          for ind in range(len(self.obs)):
               p_ind = []
               for site in range(len(self.obs[ind])):
                    p_ind_site = {}
                    for state_ind in range(len(self.states)):
                         p_ind_site[self.states[state_ind]] = \
                              eexp(gamma[ind][site][state_ind])
                    p_ind.append(p_ind_site)
               p.append(p_ind)
          return p
