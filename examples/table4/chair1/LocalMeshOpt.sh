green="\033[1;92m"
red="\033[1;91m"
yellow="\033[1;93m"
nocolor="\033[0m"
if [ "$#" -lt 3 ]; then
    echo "Usage: ./LocalMeshOpt.sh tangled.hex.<vtk> ref.<vtk|mesh|off> targetMinScaledJacobian result.<vtk> [angle(165 default)] [iters(20 default)] [smoothVolume=true]"
    echo "Example1: ./LocalMeshOpt.sh input.vtk ref.tet.vtk 0.4 result.vtk 170 20 true"
    echo "Example2: ./LocalMeshOpt.sh input.vtk input.vtk 0.5 result.vtk"
    echo "Example3: ./LocalMeshOpt.sh input.vtk input.vtk 0.6 "
    exit
fi
input=$1
ref=$2
userMSJ=$3
output="result.vtk"
if [ "$#" -gt 3 ]; then
    output=$4
fi
angle="165"
if [ "$#" -eq 5 ]; then
    angle="$5"
fi
globaliters=20
if [ "$#" -eq 6 ]; then
    angle="$5"
    globaliters=$6
fi
smoothVolume="false"
globaliters=20
if [ "$#" -eq 7 ]; then
    angle="$5"
    globaliters=$6
    smoothVolume="true"
fi

echo "input  = $input"
echo "ref    = $ref"
echo "angle  = $angle"
echo "output = $output"

mkdir LMO
cd LMO
cp ../$input ../$ref . 2>log
# This is very important
export PATH=../../../../code/bin:$PATH
echo "######## Preprocessing ########"
scaled=scaled.vtk
orig_scaled=orig_scaled.vtk
out=out.vtk

#MeshScale will create scale.txt
MeshScale $input scaled.vtk > log
if [ "$input" != "$ref" ]; then
    scale=`awk 'END {print $NF}' log`
    MeshScale $ref $orig_scaled "scale=$scale" auto=false orig=input.vtk >>log
fi
 
pi=`echo "4*a(1)" | bc -l`
deg=`echo "180 - $angle" | bc -l`
rad=`echo "$deg*($pi/180)" | bc -l`
cosangle=`echo "c($rad)" | bc -l`
#echo "cosangle=$cosangle"

if [ "$input" != "$ref" ]; then
    SmoothMap "orig=$orig_scaled" "input=$scaled" "output=$out" iters=10 preserveSharpFeature=true treatSharpFeatureAsRegular=false "cosangle=$cosangle" smoothVolume=$smoothVolume 1>>log 2>>log
else
    SmoothMap "orig=$scaled" "input=$scaled" "output=$out" iters=10 preserveSharpFeature=true treatSharpFeatureAsRegular=false "cosangle=$cosangle" 1>>log 2>>log
fi

echo "######## Untangling ########"
alpha=1000
beta=1000
xi=0.6
inverted=1000
while true
do
    echo  "alpha=$alpha beta=$beta xi=$xi"
    LocalMeshOpt  "input=$out" "alpha=$alpha" "beta=$beta" gamma=1.0 "anisotropy=$xi" allowBigStep=true minScaledJacobian=0.0 "cosangle=$cosangle" useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=$out projectToTargetSurface=true 1>>log 2>>log
    MeshQuality BestLocalOpt.vtk > quality
    inverted=`awk 'END {print $NF}' quality`
    if [ "$inverted" -ne 0 ]; then
        if (( $(echo "$xi <= 0.2" |bc -l) )); then
            let alpha="$alpha/2"
            let beta="$beta/2"
            xi=0.6
        else
            xi=`echo "$xi - 0.1" | bc -l`
        fi
    else
        break;
    fi
done

cp BestLocalOpt.vtk untangled.vtk
next=untangled.vtk
if [ "$alpha" -lt 500 ]; then
    echo "######## SLIM Deformation ########"
    SlimHex "input=$next" "target=$out" result=slim.vtk iters=20 1>>log 2>>log
    next=slim.vtk
fi

echo "######## MSJ Improving ########"
MSJ=(0.1 0.2 0.3 0.4 0.45 0.50 0.55 0.60 0.65 0.69)
#MSJ=(0.1 0.2 0.3 0.4 0.45 0.50 0.55 0.60 0.65 0.70 0.72 0.73)
#MSJ=(0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.50 0.55 0.60 0.65 0.70 0.72)
for targetMSJ in "${MSJ[@]}"
do
    if (( $(echo "$userMSJ < $targetMSJ" | bc -l) )); then
        break;
    fi 
    echo "userMSJ=$userMSJ Try targetMSJ = $targetMSJ"
    alpha=1000
    beta=1000
    xi=0.6
    obtain=0
    while [ "$obtain" -eq 0 ]
    do
        if [ "$alpha" -lt 500 ]; then
            break 2;
        fi
        echo  "alpha=$alpha beta=$beta xi=$xi"
        #echo  "LocalMeshOpt  \"input=$next\" \"alpha=$alpha\" \"beta=$beta\" gamma=1.0 \"anisotropy=$xi\" allowBigStep=true \"minScaledJacobian=$targetMSJ\" \"cosangle=$cosangle\" useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=$out projectToTargetSurface=true 1>>log 2>>log"
        LocalMeshOpt  "input=$next" "alpha=$alpha" "beta=$beta" gamma=1.0 "anisotropy=$xi" allowBigStep=true "minScaledJacobian=$targetMSJ" "cosangle=$cosangle" useSmallBlock=true blockSize=50 localIters=20 iters=$globaliters targetSurface=$out projectToTargetSurface=true 1>>log 2>>log
        MeshQuality BestLocalOpt.vtk > quality
        inverted=`awk 'END {print $NF}' quality`
        sed '3q;d' quality > msj
        msj=`awk 'END {print $NF}' msj`
        worse=`echo "$msj <= $targetMSJ" |bc -l`
        if [[ ( "$inverted" -ne 0 )  ||  ( "$worse" -ne 0 ) ]]; then
            if (( $(echo "$xi <= 0.2" |bc -l) )); then
                let alpha="$alpha/2"
                let beta="$beta/2"
                xi=0.6
            else
                xi=`echo "$xi - 0.1" | bc -l`
            fi
        else
            echo -e "${green}######## Obtained targetMSJ = $targetMSJ${nocolor}"
            cp BestLocalOpt.vtk $output
            cp BestLocalOpt.vtk "MSJ_$targetMSJ.vtk"
            next=$output
            obtain=1
        fi
    done
done

if [ -f $output ];then
    #cp $output scaled.$output
    extract_surface scaled.vtk orig.surface.off 1
    SmoothMap orig=orig.surface.off input=$output output=smooth.vtk iters=20 preserveQuality=true preserveSharpFeature=true treatSharpFeatureAsRegular=false 1>log 2>log
    MeshScale smooth.vtk mapback.vtk orig=$input
    ########################################
    MeshQuality mapback.vtk > quality
    inverted=`awk 'END {print $NF}' quality`
    sed '3q;d' quality > msj
    msj=`awk 'END {print $NF}' msj`

    MeshQuality $input > quality
    inverted=`awk 'END {print $NF}' quality`
    sed '3q;d' quality > msj
    inputmsj=`awk 'END {print $NF}' msj`
    worse=`echo "$msj <= $inputmsj" |bc -l`
    ########################################
    if [[ ( "$worse" -ne 0 ) ]]; then
        cp "../$input" "../$output"
    else
        cp mapback.vtk "../$output"
        cp mapback.vtk $output
    fi

    extract_surface $output result.surface.off 1
    extract_surface $ref ref.surface.off 1
    metro ref.surface.off result.surface.off
    echo -e "${yellow}MeshQuality $output metric=Scaled Jacobian${nocolor}"
    MeshQuality $output
    #rm BestLocalOpt.vtk Faces* Global* Local*.vtk MeshOpt*.vtk Smooth*.vtk Hex*.vtk iter*.vtk temp.vtk quality out.vtk opt.vtk msj out_Face.vtk scaled*.vtk orig_scaled.vtk 1>>log 2>>log
    #mv log ../
    cd ..
    #rm -rf LMO
else 
    ########################################
    MeshQuality $next > quality
    inverted=`awk 'END {print $NF}' quality`
    sed '3q;d' quality > msj
    msj=`awk 'END {print $NF}' msj`

    MeshQuality $input > quality
    inverted=`awk 'END {print $NF}' quality`
    sed '3q;d' quality > msj
    inputmsj=`awk 'END {print $NF}' msj`
    worse=`echo "$msj <= $inputmsj" |bc -l`
    ########################################
    if [[ ( "$worse" -ne 0 ) ]]; then
        cp "../$input" "../$output"
    else
        cp $next "../$output"
    fi
    cp $next $output
    cd ..
    rm Faces* Global* Local*.vtk MeshOpt*.vtk Smooth*.vtk quality msj out.vtk opt.vtk out_Face.vtk temp.vtk 1>>log 2>>log
    MeshQuality $output
fi


