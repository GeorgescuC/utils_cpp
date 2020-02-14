CC=g++
CFLAGS=-Wall
INCLUDE=-I/usr/local/opt/zlib/include -lhts
OPTI=-O2

all: sift_bam_max_cov sam_to_gtf
.PHONY : all

sift_bam_max_cov: sift_bam_max_cov.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(OPTI) -o sift_bam_max_cov.exe sift_bam_max_cov.cpp

sam_to_gtf: sam_to_gtf.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(OPTI) -o sam_to_gtf.exe sam_to_gtf.cpp

