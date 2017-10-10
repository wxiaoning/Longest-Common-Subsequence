## Getting started

	git clone https://github.com/wxiaoning/longest-common-subsequence.git
	cd longest-common-subsequence;
	if you're on the KNC platforms
		then use the following make command:
		make ARCH=LCS_KNC
	else if you're on the KNL platforms
		then use the following make command:
		make ARCH=LCS_KNL

## Introduction

	Longest-Common-Subsequence(LCS) is a software package for accelrating bit-parallel
longest common subsequence algorithm on Xeon Phi clusters, including KNC cluster and KNL cluster.
 The database used in the experiments are synthetic, you can use the file input.c and input_ref.c to generate synthetic subject and query data independently. You need to create a new folder named database firstly. The subject and query files are stored in database folder. The subject file is named seq.fasta and the query file is named input_seq[number].fasta. If you're running on the cluster, there're as many query files as MIC cards named input_seq0.fasta, input_seq1.fasta, input_seq2.fasta and so on. The results are stored in result[number].fasta. It gives the length of two sequences.

## Availability

	The latest source code of LCS is freely available st github. After you acquire the source code, use the front `make` instructions to compile and copy the single executable `main` to the destination you want.
