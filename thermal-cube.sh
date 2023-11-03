#!/bin/bash

# solution temperature and conductivity
T="1 + sin(2*x)^2 * cos(3*y)^2 * sin(4*z)"
k="1 + 0.5*T(x,y,z)"

# combinations
# ints="full reduced extended"
ints="full"
bcs="dirichlet neumann"
elems="tet hex"
algos="struct delaunay"
orders="1 2"
# cs="12 16 20 24 28 32 36 38 40 48 60"
cs="12 16"

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
T(x,y,z) := ${T};
k(x,y,z) := ${k};
q(x,y,z) := -(diff(k(x,y,z) * diff(T(x,y,z), x), x) + diff(k(x,y,z) * diff(T(x,y,z), y), y) + diff(k(x,y,z) * diff(T(x,y,z), z), z));
stringout("thermal-cube-q.txt", q(x,y,z));
stringout("thermal-cube-qx.txt", -k(x,y,z) * diff(T(x,y,z),x));
stringout("thermal-cube-qy.txt", -k(x,y,z) * diff(T(x,y,z),y));
stringout("thermal-cube-qz.txt", -k(x,y,z) * diff(T(x,y,z),z));
EOF


# read back the string with q(x,y) and the partial derivatives
q=$(cat thermal-cube-q.txt | tr -d ';\n')
qx=$(cat thermal-cube-qx.txt | tr -d ';\n')
qy=$(cat thermal-cube-qy.txt | tr -d ';\n')
qz=$(cat thermal-cube-qz.txt | tr -d ';\n')


# report what we found
cat << EOF
# manufactured solution
${T}
${k}
# source terms
q(x,y) = ${q}
qx(x,y) = ${qx}
qy(x,y) = ${qy}
qz(x,y) = ${qz}

EOF


# empty ppls
rm -f thermal-cube-fits.ppl
echo "plot \\" > thermal-cube-einf.ppl
echo "plot \\" > thermal-cube-e2.ppl

for bc in ${bcs}; do
 for int in ${ints}; do
 
  echo "plot \\" > thermal-cube-${bc}-${int}-einf.ppl
  echo "plot \\" > thermal-cube-${bc}-${int}-e2.ppl

  for elem in ${elems}; do
   for order in ${orders}; do
    for algo in ${algos}; do

     # do not compute hex-2
     #                reduced-tet-2
#      if [[ ( "x${elem}" != "xhex" || "x${order}" != "x2" ) &&
#            ( "x${elem}" != "xtet" || "x${order}" != "x2" || "x${int}" != "xreduced" )  ]]; then
      dat="thermal-cube-${bc}-${int}-${elem}-${order}-${algo}"
      rm -f ${dat}.dat
      echo ${dat}
      echo "-----------------------------------------------------------"
     
      for c in ${cs}; do
       name="thermal-cube-${bc}-${int}-${elem}-${order}-${algo}-${c}"
      
       if [ ! -e cube-${elem}-${algo}-${c}.msh ]; then
        if [[ "x${elem}" = "xhex" && "x${algo}" = "xdelaunay" ]]; then
          lc=$(echo "PRINT 1/(1+(${c}-12)*24/(64-12))" | feenox -)
        else  
          lc=$(echo "PRINT 1/${c}" | feenox -)
        fi
        gmsh -v 0 -3 cube.geo ${elem}.geo ${algo}.geo -clscale ${lc} -o cube-${elem}-${algo}-${c}.msh
       fi
       
       # prepare input JSON
       if [ ! -e ${name}.json ]; then
        cp thermal-cube_${bc}.json ${name}.json
        sed -i "s/\$m\\$/cube-${elem}-${algo}-${c}.msh/g" ${name}.json
        sed -i "s/\$i\\$/${int^}/g" ${name}.json
        sed -i "s/\$o\\$/${order}/g" ${name}.json
        sed -i "s/\$T\\$/${T}/g" ${name}.json
        sed -i "s/\$qx\\$/${qx}/g" ${name}.json
        sed -i "s/\$qy\\$/${qy}/g" ${name}.json
        sed -i "s/\$qz\\$/${qz}/g" ${name}.json
        sed -i "s!\$q!${q}!" ${name}.json
        
        # reflex uses just T not T(x,y,z(
        kTnoargs=$(echo ${k} | sed 's/T(x,y,z)/T/g')
        sed -i "s/\$k\\$/${kTnoargs}/g" ${name}.json
       fi
       
       # call reflex (if needed, otherwise we use a cached result)
       if [ ! -e output/${name}_result_1.vtu ]; then
        reflexCLI -i ${name}.json > /dev/null
       fi
     
       # compute errors out of the vtu
       if [ -e output/${name}_result_1.vtu ]; then
        # o warns about different size in python and in C
        meshio convert output/${name}_result_1.vtu --output-format=gmsh ${name}.msh 2> /dev/null
        feenox error3d.fee ${name} "${T}" ${c} | tee -a ${dat}.dat
       fi
       
      done  # c 

      feenox fit.fee ${dat} | sed 's/-/_/g' | sed 's/e_0/e-0/' | sed 's/* _/* -/' >> thermal-cube-fits.ppl
     
      cat << EOF >> thermal-cube-einf.ppl
     "${dat}.dat"                              u 1:2 w lp pt ${pt[${elem}${order}]} lw ${lw[${bc}]} lt ${lt[${algo}]} color ${co[${int}]}  ti "${bc}-${int}-${elem}${order}-${algo} = " + e_inf_thermal_cube_${bc}_${int}_${elem}_${order}_${algo}_title,\\
EOF
      cat << EOF >> thermal-cube-e2.ppl
     "${dat}.dat"                              u 1:3 w lp pt ${pt[${elem}${order}]} lw ${lw[${bc}]} lt ${lt[${algo}]} color ${co[${int}]}  ti "${bc}-${int}-${elem}${order}-${algo} = " + e_2_thermal_cube_${bc}_${int}_${elem}_${order}_${algo}_title,\\
EOF
     

      cat << EOF >> thermal-cube-${bc}-${int}-einf.ppl
     "${dat}.dat"                              u 1:2 w lp pt ${pt[${elem}${order}]} lw ${lw[${bc}]} lt ${lt[${algo}]} color ${co[${int}]}  ti "${bc}-${int}-${elem}${order}-${algo} = " + e_inf_thermal_cube_${bc}_${int}_${elem}_${order}_${algo}_title,\\
EOF
      cat << EOF >> thermal-cube-${bc}-${int}-e2.ppl
     "${dat}.dat"                              u 1:3 w lp pt ${pt[${elem}${order}]} lw ${lw[${bc}]} lt ${lt[${algo}]} color ${co[${int}]}  ti "${bc}-${int}-${elem}${order}-${algo} = " + e_2_thermal_cube_${bc}_${int}_${elem}_${order}_${algo}_title,\\
EOF
     
#      fi    # not unstructured hexes
    done   # algo
   done    # order
  done     # elem
 done      # int
done       # bc

cat << EOF >> thermal-cube-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF



cat << EOF >> thermal-cube-dirichlet-full-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-dirichlet-reduced-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-dirichlet-extended-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-neumann-full-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-neumann-reduced-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-neumann-extended-einf.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF



cat << EOF >> thermal-cube-dirichlet-full-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-dirichlet-reduced-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-dirichlet-extended-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-neumann-full-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-neumann-reduced-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF

cat << EOF >> thermal-cube-neumann-extended-e2.ppl
 x**2    w l lt 3 lw 4 color gray ti "\$h^2\$",\\
 x**3    w l lt 4 lw 4 color gray ti "\$h^3\$"
EOF
