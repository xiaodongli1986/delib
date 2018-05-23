
#DElibPATH=/home/xue/workspace/delib
DElibPATH=/home/xiaodongli/software/delib

# libarary
export LIBRARY_PATH=$DElibPATH/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$DElibPATH/lib:$LD_LIBRARY_PATH

# moduls path
export demods=$DElibPATH/mods

# link libarary, include modules
export delm=-lde\ -I$demods
