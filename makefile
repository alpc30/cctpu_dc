#
# Variable of CSPPLAN
#
SRC = src/

#change
all:
	@echo "Building....optimise-cctpu"
	@cd $(SRC); make all
	@echo "Build successfull..."

clean:
	@echo "Cleaning optimise-cctpu"
	@cd $(SRC); make clean



