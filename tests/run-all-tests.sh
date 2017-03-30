# Run all tests

m=../meshgeometry_mac

echo -n "Test 1. Translate: "
# run the test
$m -i input/cube.ply -translate 10 10 10 -o output_computed/translate-test_cube.ply
# check result
a=$($m -i output_correct/translate-test_cube.ply -printBarycentre)
b=$($m -i output_computed/translate-test_cube.ply -printBarycentre)
if [ "$a" == "$b" ]; then
	echo "Passed."
else
	echo "ERROR."
fi

