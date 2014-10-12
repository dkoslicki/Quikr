# ==============================================================================
# Quikr.jl
#
# Authors: David Koslicki (david.koslicki@math.oregonstate.edu)
#
# Takes a fasta file as input and performs the Quikr Algorithm on it.
# ==============================================================================

using ArgParse
using HDF5
include("lsqnonneg.jl")
include("ConvertToCAMIOutput.jl")


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
		"--kmer_count_path", "-k"
			help = "Full path to the kmer_total_count file obtained from https://github.com/mutantturkey/dna-utils"
			default = "kmer_total_count"
		"--training_file", "-t"
			help = "Training database to use. Choices are: SEK, Quikr. SEK is more accurate, but slower, whereas Quikr is faster, but less accurate. Default is Quikr"
			default = "Quikr"
			
    end
    return parse_args(s)
end


#Parse the args
parsed_args = parse_commandline()
input_file = parsed_args["input_file"]
lambda = int(parsed_args["lambda"])
output_file = parsed_args["output_file"]
kmer_count_path = parsed_args["kmer_count_path"]
training_file = parsed_args["training_file"]

#Form the 6mer counts
counts = map(x->int(strip(x)),readlines(`$kmer_count_path -i $input_file -k 6 -c`));

#normalize the counts
counts = counts/sum(counts);

if training_file == "Quikr"
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
	output_level = 0; #Since we don't have hypothetical organisms
	ConvertToCAMIOutput(x, "../../data/trainset7_taxonomy.txt", output_level, output_file)
elseif training_file == "SEK"
	#Read in the training database
	A = h5read("../../data/trainset7_112011_allseqslongerthan700-SEKTrainingMatrix-bitShift100-windowLength400-N6C.h5","/data");
	
	#Read in the block matrix to return to the trainset7_SEK basis
	blockMatrix = h5read("../../data/trainset7_112011_allseqslongerthan700-SEKTrainingMatrix-bitShift100-windowLength400-blockMatrix.h5","/data");

	#Form the Aaux
	Aaux = [ones(1,size(A,2)); lambda*A];
	yaux = [0;lambda*counts];

	#Perform the reconstruction
	x = lsqnonneg(Aaux, yaux);
	
	#Return to the trainset7_SEK basis
	x = blockMatrix*x;

	#Normalize the output
	x = x/sum(x);

	#Write the output to file
	output_level = 0; #Since we don't have hypothetical organisms
	ConvertToCAMIOutput(x, "../../data/trainset7_SEK_taxonomy.txt", output_level, output_file)
end










#How I saved the training data
#fid = h5open("trainset7_112011N6C.h5","w");
#fid["/data", "chunk", (512,10046), "compress", 7] = Atemp;
#close(fid)

