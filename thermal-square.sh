#!/bin/bash -e

# proposed solution temperature and conductivity
T="1 + sin(2*x)^2 * cos(3*y)^2"
k="1 + x - 0.5*y"

# combinations
ints="full reduced extended"
bcs="dirichlet neumann"
elems="tri quad"
algos="struct frontal delaunay"
orders="1 2"
# cs="8 10 12 14 16 20 24 28 32 36 40 48 52 56 58 60"
cs="8 10 12 14 16"
# --------------------------------------

# sanity check
for i in reflexCLI feenox gmsh maxima meshio tr sed; do
 if [ -z "$(which ${i})" ]; then
  echo "error: ${i} not installed"
  exit 1
 fi
done


. styles.sh


# call maxima to compute source term
maxima --very-quiet << EOF > /dev/null
T(x,y) := ${T};
k(x,y) := ${k};
q(x,y) := -(diff(k(x,y) * diff(T(x,y), x), x) + diff(k(x,y) * diff(T(x,y), y), y));
stringout("thermal-square-q.txt", q(x,y));
stringout("thermal-square-qx.txt", -k(x,y) * diff(T(x,y),x));
stringout("thermal-square-qy.txt", -k(x,y) * diff(T(x,y),y));
EOF



# read back the string with q(x,y) and the partial derivatives
q=$(cat thermal-square-q.txt | tr -d ';\n')
qx=$(cat thermal-square-qx.txt | tr -d ';\n')
qy=$(cat thermal-square-qy.txt | tr -d ';\n')


# report what we found
cat << EOF
# manufactured solution
${T}
${k}
# source terms
q(x,y) = ${q}
qx(x,y) = ${qx}
qy(x,y) = ${qy}

EOF


# empty ppls
rm -f thermal-square-fits.ppl
echo "plot \\" > thermal-square-einf.ppl
echo "plot \\" > thermal-square-e2.ppl

for bc in ${bcs}; do
 for int in ${ints}; do
 
  echo "plot \\" > thermal-square-${bc}-${int}-einf.ppl
  echo "plot \\" > thermal-square-${bc}-${int}-e2.ppl

  for elem in ${elems}; do
   for order in ${orders}; do
    for algo in ${algos}; do

     if [[ "x${elem}" != "xtri" || "x${int}" != "xreduced" || "x${order}" != "x2" ]]; then
    
     dat="thermal-square-${bc}-${int}-${elem}-${order}-${algo}"
     rm -f ${dat}.dat
     echo ${dat}
     echo "-----------------------------------------------------------"
     
     for c in ${cs}; do
      name="thermal-square-${bc}-${int}-${elem}-${order}-${algo}-${c}"
      
      if [ ! -e square-${elem}-${algo}-${c}.msh ]; then
       lc=$(echo "PRINT 1/${c}" | feenox -)
       gmsh -v 0 -2 square.geo ${elem}.geo ${algo}.geo -clscale ${lc} -o square-${elem}-${algo}-${c}.msh
      fi
      
      # prepare input JSON
      if [ ! -e ${name}.json ]; then
       cp thermal-square_${bc}.json ${name}.json
       sed -i "s/\$m\\$/square-${elem}-${algo}-${c}.msh/g" ${name}.json
       sed -i "s/\$i\\$/${int^}/g" ${name}.json
       sed -i "s/\$o\\$/${order}/g" ${name}.json
       sed -i "s/\$T\\$/${T}/g" ${name}.json
       sed -i "s/\$qx\\$/${qx}/g" ${name}.json
       sed -i "s/\$qy\\$/${qy}/g" ${name}.json
       sed -i "s/\$q\\$/${q}/g" ${name}.json
       sed -i "s/\$k\\$/${k}/g" ${name}.json
      fi
      
      # call reflex (if needed, otherwise we use a cached result)
      if [ ! -e output/${name}_result1.vtu ]; then
       reflexCLI -i ${name}.json > /dev/null
      fi
    
      # compute errors out of the vtu
      if [ -e output/${name}_result_1.vtu ]; then
       # meshio warns about different size in python and in C
       meshio convert output/${name}_result_1.vtu --output-format=gmsh ${name}.msh 2> /dev/null
       feenox error2d.fee ${name} "${T}" ${c} | tee -a ${dat}.dat
      fi  
      
     done  # c 

     feenox fit.fee ${dat} | sed 's/-/_/g' | sed 's/e_0/e-0/' | sed 's/* _/* -/' >> thermal-square-fits.ppl
     
    cat << EOF >> thermal-square-einf.ppl
     "${dat}.dat"                              u 1:2 w lp pt ${pt[${elem}${order}]} lw ${lw[${bc}]} lt ${lt[${algo}]} color ${co[${int}]}  ti "${bc}-${int}-${elem}${order}-${algo} = " + e_inf_thermal_square_${bc}_${int}_${elem}_${order}_${algo}_title,\\
EOF
    cat << EOF >> thermal-square-e2.ppl
     "${dat}.dat"                              u 1:3 w lp pt ${pt[${elem}${order}]} lw ${lw[${bc}]} lt ${lt[${algo}]} color ${co[${int}]}  ti "${bc}-${int}-${elem}${order}-${algo} = " + e_2_thermal_square_${bc}_${int}_${elem}_${order}_${algo}_title,\\
EOF
     

  cat << EOF >> thermal-square-${bc}-${int}-einf.ppl
     "${dat}.dat"                              u 1:2 w lp pt ${pt[${elem}${order}]} lw ${lw[${bc}]} lt ${lt[${algo}]} color ${co[${int}]}  ti "${bc}-${int}-${elem}${order}-${algo} = " + e_inf_thermal_square_${bc}_${int}_${elem}_${order}_${algo}_title,\\
EOF
  cat << EOF >> thermal-square-${bc}-${int}-e2.ppl
     "${dat}.dat"                              u 1:3 w lp pt ${pt[${elem}${order}]} lw ${lw[${bc}]} lt ${lt[${algo}]} color ${co[${int}]}  ti "${bc}-${int}-${elem}${order}-${algo} = " + e_2_thermal_square_${bc}_${int}_${elem}_${order}_${algo}_title,\\
EOF
     
     fi
    done   # algo
   done    # order
  done     # elem
 done      # int
done       # bc

cat << EOF >> thermal-square-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-square-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF



cat << EOF >> thermal-square-dirichlet-full-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-square-dirichlet-reduced-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

# cat << EOF >> thermal-square-dirichlet-extended-einf.ppl
#  x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
#  x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
# EOF

cat << EOF >> thermal-square-neumann-full-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-square-neumann-reduced-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

# cat << EOF >> thermal-square-neumann-extended-einf.ppl
#  x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
#  x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
# EOF



cat << EOF >> thermal-square-dirichlet-full-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-square-dirichlet-reduced-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

# cat << EOF >> thermal-square-dirichlet-extended-e2.ppl
#  x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
#  x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
# EOF

cat << EOF >> thermal-square-neumann-full-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-square-neumann-reduced-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

# cat << EOF >> thermal-square-neumann-extended-e2.ppl
#  x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
#  x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
# EOF
