# If not running interactively, don't do anything
case $- in
    *i*) ;;
      *) return;;
esac

# Get setup for BOUT++
source /work/e281/e281/$USER/BOUT-configs/archer/bout.env
