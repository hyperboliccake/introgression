import sys
import sim_process

# takes in two sets of introgression calls on same set of data and
# compares them (for example, actual and predicted regions, or calls
# from two different prediction methods)

fn1 = sys.argv[1]
fn2 = sys.argv[2]

f1 = open(fn1, 'r')
f2 = open(fn2, 'r')

line1 = f1.readline()
line2 = f2.readline()

while line1 != '':

    d1, rep1, line1 = sim_process.read_introgression(f1, line1)
    d2, rep2, line2 = sim_process.read_introgression(f2, line2)
    assert rep1 == rep2, str(rep1) + ' ' + str(rep2)

    

assert line2 == '', line2
