docker run -d -v $1/data:/data -v $1/results:/results sequencegenie \
pathway \
/results \
/data \
https://ice.synbiochem.co.uk \
$2 \
$3 \
gaattcaaaagatcttttaagaag \
ttactcgagtttggatcc \
True