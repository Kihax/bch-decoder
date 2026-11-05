PLATFORM = /home/killian/Xilinx/2025.2/Vitis/base_platforms/vek385_base/vek385_base.xpfm
WORK_DIR = Work

all: compile simulate

compile:
	aiecompiler --target=hw \
	            --platform=$(PLATFORM) \
	            --include="." \
	            --workdir=$(WORK_DIR) \
	            main.cpp

simulate:
	aiesimulator --pkg-dir=$(WORK_DIR)

clean:
	rm -rf $(WORK_DIR) aiesimulator_output .AIE_SIM_CMD_history *.log vitis_*