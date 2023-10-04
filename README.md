# masked_bnn_riscv
Efficiently masked binarised neural network for RISC-V embedded devices

We use SiFive freedom studio to create our execution environment.
To run our code, create a RISC-V project for the HiFive1 Rev B in Freedom Studio, and replace the src folder with the one in this repo.
The run/debug tooling in Freedom Studio will link/compile and flash the binary onto the device.

It may be needed to update the linker script to allocate a larger stack and heap section in order to make a binary that fits the constraints of the device.
We include the linker script that works for our device and source code.
