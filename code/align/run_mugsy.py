import os

drc = '/net/akey/vol2/aclark4/nobackup/100_genomes/genomes/'
refc = 'S288c_SGD-R64.fa'

drp = '/net/akey/vol2/aclark4/nobackup/100_genomes/paradoxus/strains/CBS432/assembly/'
refp = 'genome.fa'

dc = '/net/akey/vol2/aclark4/nobackup/100_genomes/genomes/'
c = ['Sigma1278b.fa', 'SK1.fa', 'yjm1078.fa', 'yjm1083.fa', 'yjm1129.fa', 'yjm1133.fa', 'yjm1190.fa', 'yjm1199.fa', 'yjm1202.fa', 'yjm1208.fa', 'yjm1242.fa', 'yjm1244.fa', 'yjm1248.fa', 'yjm1250.fa', 'yjm1252.fa', 'yjm1273.fa', 'yjm1304.fa', 'yjm1307.fa', 'yjm1311.fa', 'yjm1326.fa', 'yjm1332.fa', 'yjm1336.fa', 'yjm1338.fa', 'yjm1341.fa', 'yjm1342.fa', 'yjm1355.fa', 'yjm1356.fa', 'yjm1381.fa', 'yjm1383.fa', 'yjm1385.fa', 'yjm1386.fa', 'yjm1387.fa', 'yjm1388.fa', 'yjm1389.fa', 'yjm1399.fa', 'yjm1400.fa', 'yjm1401.fa', 'yjm1402.fa', 'yjm1415.fa', 'yjm1417.fa', 'yjm1418.fa', 'yjm1419.fa', 'yjm1433.fa', 'yjm1434.fa', 'yjm1439.fa', 'yjm1443.fa', 'yjm1444.fa', 'yjm1447.fa', 'yjm1450.fa', 'yjm1460.fa', 'yjm1463.fa', 'yjm1477.fa', 'yjm1478.fa', 'yjm1479.fa', 'yjm1526.fa', 'yjm1527.fa', 'yjm1549.fa', 'yjm1573.fa', 'yjm1574.fa', 'yjm1592.fa', 'yjm1615.fa', 'yjm189.fa', 'yjm193.fa', 'yjm195.fa', 'yjm244.fa', 'yjm248.fa', 'yjm270.fa', 'yjm271.fa', 'yjm320.fa', 'yjm326.fa', 'yjm428.fa', 'yjm450.fa', 'yjm451.fa', 'yjm453.fa', 'yjm456.fa', 'yjm470.fa', 'yjm541.fa', 'yjm554.fa', 'yjm555.fa', 'yjm627.fa', 'yjm681.fa', 'yjm682.fa', 'yjm683.fa', 'yjm689.fa', 'yjm693.fa', 'yjm969.fa', 'yjm972.fa', 'yjm975.fa', 'yjm978.fa', 'yjm981.fa', 'yjm984.fa', 'yjm987.fa', 'yjm990.fa', 'yjm993.fa', 'yjm996.fa']

cmd_string = ''

cmd_string += 'export MUGSY_INSTALL=~/software/mugsy; '
cmd_string += 'export PATH=$PATH:$MUGSY_INSTALL:$MUGSY_INSTALL/mapping; '
cmd_string += 'export PERL5LIB=$MUGSY_INSTALL/perllibs; '

chrms_roman = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XIV']

for x in c:

    cmd_string = ''
    
    cmd_string += 'export MUGSY_INSTALL=~/software/mugsy; '
    cmd_string += 'export PATH=$PATH:$MUGSY_INSTALL:$MUGSY_INSTALL/mapping; '
    cmd_string += 'export PERL5LIB=$MUGSY_INSTALL/perllibs; '
    
    for chrm in chrms_roman:
        cmd_string += '~/software/mugsy/mugsy --directory ~/introgression/100genomes/chrm/alignments/ --prefix S288c_CBS432_' + x[:-3] + '_chr' + chrm + ' ' + drc + refc[:-3] + '/' + refc[:-3] + '_chr' + chrm + '.fa' + ' ' + drp + 'CBS432/CBS432_chr' + chrm + '.fa' + ' ' + dc + x[:-3] + '/' + x[:-3] + '_chr' + chrm + '.fa' + '; '
    print x
    # commands can only be up to a certain length so break it up this way
    os.system(cmd_string)
