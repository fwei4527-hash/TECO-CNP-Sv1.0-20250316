#!/bin/bash

# -----------------------------------------------------------------------
# Compile and run the TECO-CNP biogeochemical model
# You need to install the dependences according to the README file
# -----------------------------------------------------------------------

# Step 1: Compile the Fortran source files
echo "Compiling Fortran source files..."
gfortran -ffree-form -c FileSize.f90 ParasModule.f90 SASpinUp.f90 LIMITATION.f90 NPUptakeDemand.f90 NPDynamic.f90 MCMC.f90 TECO_CNP_main.f90

# Step 2: Link the object files into an executable
echo "Linking object files..."
#gfortran -o teco_cnp.exe FileSize.o ParasModule.o SASpinUp.o LIMITATION.o NPUptakeDemand.o NPDynamic.o MCMC.o TECO_CNP_main.o -L/YourPath -llapacke -llapack -lrefblas
#gfortran -o teco_cnp.exe FileSize.o ParasModule.o SASpinUp.o LIMITATION.o NPUptakeDemand.o NPDynamic.o MCMC.o TECO_CNP_main.o -L/home/wan_fx/usr/lib -llapacke -llapack -lblas
gfortran -o teco_cnp.exe FileSize.o ParasModule.o SASpinUp.o LIMITATION.o NPUptakeDemand.o NPDynamic.o MCMC.o TECO_CNP_main.o -llapacke -llapack -lblas

# Step 3: Execute the program with the provided arguments
# To run the TECO-CNP model, you can adjust the configuration by using different combinations of arguments. 
# Please refer to the definitions of the arguments in the `TECO-CNP_main.f90` file.
# Validation check: SPINUP must be either 1 or nspinup should be >= 1
# WARNING: These functions cannot be executed concurrently
# IMPORTANT: Please verify file permissions before proceeding
# 参数含义：
# 1: start_year  = 2001
# 2: end_year    = 2010
# 3: CYCLE_CNP   = 3     (CNP 模式)
# 4: MCMC        = 0     (不用 MCMC)
# 5: NDSPINUP    = 1     (开启 spinup)
# 6: nspinup     = 500   (spinup 循环次数)
# 7: SensTest    = 0     (不用敏感性分析)

# echo "Running the executable..."
# #./teco_cnp.exe 2022 2024 3 0 0 5000 0
# ./teco_cnp.exe 2021 2024 3 0 1 2000 0 0 0       


# echo "Done."

 ----------------------------------------------
echo "Running P addition treatments..."

for PADD in 0 25 50 100
do
    echo "======================================"
    echo "Running P addition: ${PADD} kg P ha-1 yr-1"
    echo "======================================"

    rm -rf ../output/sim/teco_cnp
    mkdir -p ../output/sim/teco_cnp

    ./teco_cnp.exe 2021 2024 3 0 0 5000 0 1 ${PADD}

    if [ ! -d ../output/sim/teco_cnp ]; then
        echo "ERROR: ../output/sim/teco_cnp was not created."
        echo "Model likely failed for PADD=${PADD}"
        exit 1
    fi

    mkdir -p "../output/sim/${PADD}磷添加"
    rm -rf "../output/sim/${PADD}磷添加/teco_cnp"
    mv ../output/sim/teco_cnp "../output/sim/${PADD}磷添加/teco_cnp"

done

echo "Done."
