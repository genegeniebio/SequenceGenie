docker run -d -v $1/data:/data -v $1/results:/results sequencegenie \
pathway \
/results \
"/data/barcodes.csv" \
"/data" \
"/data/ice_ids.txt" \
https://ice.synbiochem.co.uk \
neil.swainston@manchester.ac.uk \
$2 \
gaattcaaaagatcttttaagaagg \
ttactcgagtttggatcc \
True