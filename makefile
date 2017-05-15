# No fancy makefile here -- just keeping it simple.
# You'll need to download gf_complete, compile it, 
# and then copy over gf_complete.h and gf_complete.a

all: sd_code generate_random_data_file erase_data
clean: 
	rm -f *.o sd_code generate_random_data_file erase_data

sd_code: sd_code.c gf_complete.h gf_complete.a
	gcc -O3 -o sd_code sd_code.c gf_complete.a

generate_random_data_file: generate_random_data_file.c
	gcc -O3 -o generate_random_data_file generate_random_data_file.c

erase_data: erase_data.c
	gcc -O3 -o erase_data erase_data.c

