if [ "$#" -ne 5 ]; then
    echo "Usage: $0 molarmass (in g/mol) density (kg/m3) length lx (m) ly (m) lz (m)" >&2
    echo " "
    echo "Example: $0 18.01534 1000 10e-10 10e-10 10-e10"
    echo "The above example will calculate number of molecules that \
could fit a cubic volume of length 10e-10 m for water with \
molar mass of 18.01534 g/mol and density of 1000 kg/m3"
exit 1
fi

NA=6.0221415e23
mass=$1
density=$2
Lx=$3
Ly=$4
Lz=$5
echo "Molar Mass (g/mol) = $mass"
echo "Density (kg/m^3)   = $density"
echo "Lx x Ly x Lz       = $Lx $Ly $Lz"
echo "N                  = 6.022e23*1e3*(Density/MolarMass)*Lx*Ly*Lz"
awk -v NA=$NA -v m=$mass -v d=$density -v Lx=$Lx -v Ly=$Ly -v Lz=$Lz  'BEGIN { printf "No of molecules    = %g\n", NA*(d*1e3/m)*(Lx*Ly*Lz) }'
