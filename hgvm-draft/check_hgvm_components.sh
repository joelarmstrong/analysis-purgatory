# Checks a HAL file produced with the "cycle-free isolated components"
# functionality of Cactus. Checks if the non-alt chrs are indeed isolated.
set -o errexit
set -o nounset
set -o pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 halFile genome"
    exit 1
fi

HAL_FILE=$1
GENOME=$2
MYTMP=$(mktemp -d)
# halStats displays sequences separated by ,
IFS=,
for chr in $(halStats --sequences $GENOME $HAL_FILE); do
    if echo $chr | egrep '_alt$' >/dev/null; then
        echo "Skipping alt $chr"
        continue
    fi
    echo "Testing chr $chr"
    IN_FILE=$MYTMP/in.bed
    halStats --bedSequences $GENOME $HAL_FILE | egrep "^$chr\s" > $IN_FILE
    # Lift up
    PARENT=$(halStats --parent $GENOME $HAL_FILE)
    PARENT_FILE=$MYTMP/parent.bed
    halLiftover $HAL_FILE $GENOME $IN_FILE $PARENT $PARENT_FILE
    # Lift back down
    OUT_FILE=$MYTMP/out.bed
    halLiftover $HAL_FILE $PARENT $PARENT_FILE $GENOME $OUT_FILE
    # Check output
    #
    # uncomment to show an unexpected input will actually trigger failure
    #echo "chr22_GL383582v2 3030 3101022030" >>$OUT_FILE
    if egrep -v "^(\w+_alt|$chr)\s" $OUT_FILE >$MYTMP/failure.out; then
        echo "FAILURE on initial lifting from $chr."
        cat $MYTMP/failure.out
        exit 1
    fi
    for mappedChr in $(cut -f 1 $OUT_FILE | egrep -v "^$chr$" | sort | uniq | tee test.txt | tr '\n' ','); do
        echo "Checking mapped chr $mappedChr"
        IN_FILE=$MYTMP/in.bed
        halStats --bedSequences $GENOME $HAL_FILE | egrep "^$mappedChr\s" > $IN_FILE
        # Lift up
        PARENT=$(halStats --parent $GENOME $HAL_FILE)
        PARENT_FILE=$MYTMP/parent.bed
        halLiftover $HAL_FILE $GENOME $IN_FILE $PARENT $PARENT_FILE
        # Lift back down
        OUT_FILE=$MYTMP/out.bed
        halLiftover $HAL_FILE $PARENT $PARENT_FILE $GENOME $OUT_FILE
        if egrep -v "^(\w+_alt|$chr)\s" $OUT_FILE >$MYTMP/failure.out; then
            echo "FAILURE on a secondary lifting from mapped chr $mappedChr."
            cat $MYTMP/failure.out
            exit 1
        fi
    done
done
rm -fr $MYTMP
