# ==============================================================================
# Reconstruction.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Takes a fasta file as input and performs the Quikr Algorithm on it.
# ==============================================================================

using ArgParse
using HDF5
include("lsqnonneg.jl")

#Parse arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input_file", "-i"
			help = "Fasta file to perform the reconstruction on"
        "--lambda", "-l"
			help = "Lambda parameter. The lower it is, the sparser the reconstruction. The higher it is, the more closely the kmer counts will be fit. Default is 10,000"
			default = 10000
		"--output_file", "-o"
			help = "Output text file"
    end
    return parse_args(s)
end


#Parse the args
parsed_args = parse_commandline()
input_file = parsed_args["input_file"]
lambda = parsed_args["lambda"]
output_file = parsed_args["output_file"]

#Form the 6mer counts
counts = map(x->int(strip(x)),readlines(`kmer_total_count -i $input_file -k 6 -c`));

#normalize the counts
counts = counts/sum(counts);

#Read in the training database
A = h5read("../../data/trainset7_112011N6C.h5","/data");

#Form the Aaux
Aaux = [ones(1,size(A,2)); lambda*A];
yaux = [0;lambda*counts];

#Perform the reconstruction
x = lsqnonneg(Aaux, yaux);

#Normalize the output
x = x/sum(x);

#Write the output to file
fid = open(output_file,"w")
for i = 1:length(x)
	write(fid,"$(x[i])\n")
end










#fid = h5open("trainset7_112011N6C.h5","w");
#fid["/data", "chunk", (512,10046), "compress", 7] = Atemp;
#close(fid)

