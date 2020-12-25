# WM1D: Solver for 1D Wannier-Mott model

WM1D is an open-source software to describe optical responses of excitonic systems based on the one-dimensional Wanner-Mott model.

## License

WM1D is available under GPL License.


    Copyright 2020 Shunsuke A. Sato

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Install and Run

   One can simply install WM1D by typing the following in Linux command-line:  

```
gfortran main.f90 -o WM1D -llapack
```

After the compilation, the executable file ```WM1D``` is generated.
One can run ```WM1D``` by typing the following in Linux command-line:  
```
./WM1D <input
```
Here ```input``` is the input file for the applied laser parameters. An example input file is as follows:

```
1d-2   1d-4  ! A0_pump, E0_probe
16d0   1d0   ! Tpump_fs, Tprobe_fs
1.55d0  55d0   !omega_pump_ev, omega_probe_ev
0d0          ! Tdelay_fs
```

The first column is used for amplitudes of pump and probe laser fields in the atomic unit. The component is the amplitude of the vector potential of the pump field, while the second component is the amplitude of the electric field of the probe field.  
In the second column, the first and second components are the pulse durations of the pump and probe fields, respectively. Here the values should be given in femtosecond.  
In the third column, the first and second components are the mean photon energies of the pump and probe fields, respectively. Here the values should be given in electron-volt.  
The fourth column is the time-delay between the pump and probe pulses in the unit of femtosecond.  
All the other parameters are hard-coded in ```main.f90```. For the details, see the source file ```main.f90```.