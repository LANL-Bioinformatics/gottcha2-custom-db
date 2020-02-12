GOTTCHA_DB_OBJS = options.o \
	mpi_util.o update.o parse_mapping_file.o io.o \
	fasta.o genbank.o memory_util.o write_targets.o \
	digest.o taxa.o subtract_sequence.o
	
CC = mpic++

PROFILE = #-pg
OPENMP = -fopenmp
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++0x

INC = -I. -I/home/jgans/zlib/include
LIBS = -lm /home/jgans/zlib/lib/libz.a

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(INC) -c $<

all: gottcha_db
	
gottcha_db : $(GOTTCHA_DB_OBJS) gottcha_db.o
	$(CC) $(PROFILE) -o gottcha_db $(GOTTCHA_DB_OBJS) gottcha_db.o $(LIBS) $(OPENMP)
	
clean:
	-rm -f *.o


