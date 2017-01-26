#input_filename="input/test.align"
input_filename="input/valignments.fa"
output_filename="output/test_model.csv"

python2 -m unittest unittests
./main.py -i ${input_filename} -o ${output_filename} -mf k_neighbour
# ./main.py -i ${input_filename} -o ${output_filename} -mf trivial 
