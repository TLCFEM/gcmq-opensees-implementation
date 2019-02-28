model BasicBuilder -ndm 2 -ndf 3
node        1      0 0
node        2      1 0
node        3      1 1
node        4      0 1

nDMaterial ElasticIsotropic 2 1000 .2
nDMaterial PlaneStress 1 2

element GCMQ 1 1 2 3 4 .1 1 1

fix 1 1 1 1
fix 2 1 1 1

pattern Plain 101 Constant {
load 3 0.0 1 0.0
}

constraints Plain
numberer Plain
system FullGeneral
test NormDispIncr 1E-8 10
algorithm Newton
analysis Static

analyze 1

print node 3

exit