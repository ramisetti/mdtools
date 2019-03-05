if [ "$#" -ne 5 ] && [ "$#" -ne 3 ]; then
    echo "Usage 1: $0 molarmass (in g/mol) density (kg/m3) length lx (m) ly (m) lz (m)" >&2
    echo "      2: $0 molarmass (in g/mol) density (kg/m3) Number_molecules (-)" >&2
    echo " "
    echo "Example 1: $0 18.01534 1000 10e-10 10e-10 10e-10"
    echo "The above example will calculate number of molecules that \
could fit a cubic volume of length 10e-10 m for water with \
molar mass of 18.01534 g/mol and density of 1000 kg/m3"
    echo "Example 2: $0 18.01534 1000 1500"
    echo "The above example will calculate length of a cubic box \
for 1500 water molecules with molar mass of 18.01534 g/mol and density of 1000 kg/m3"

exit 1
fi

NA=6.0221415e23
mass=$1
density=$2

if [ "$#" == 5 ]; then
    Lx=$3
    Ly=$4
    Lz=$5
    echo "Molar Mass (g/mol) = $mass"
    echo "Density (kg/m^3)   = $density"
    echo "Lx x Ly x Lz       = $Lx $Ly $Lz"
    echo "N                  = 6.0221415e23*1e3*(Density/MolarMass)*Lx*Ly*Lz"
    awk -v NA=$NA -v m=$mass -v d=$density -v Lx=$Lx -v Ly=$Ly -v Lz=$Lz  'BEGIN { printf "No of molecules    = %6.9e\n", NA*(d*1e3/m)*(Lx*Ly*Lz) }'
fi

if [ "$#" == 3 ]; then
    N=$3
    echo "Molar Mass (g/mol) = $mass"
    echo "Density (kg/m^3)   = $density"
    echo "Lx x Ly x Lz       = $Lx $Ly $Lz"
    echo "Lx                  = ((N*MolarMass/(6.0221415e23*Density)))^(1/3)"
    awk -v NA=$NA -v m=$mass -v d=$density -v n=$N 'BEGIN { printf "Length of cubic box in Angstrom   = %6.9e\n",  (((n*m)/(d*1e3*NA))^(1/3)) }'
fi
