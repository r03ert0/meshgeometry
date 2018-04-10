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

echo -n "Test 2. Edge length mean and SD of icosahedron.ply: "
# run the test
a=$($m -i input/icosahedron.ply -edgeLength)
# check result
if [ "$a" == "edgeLength: 1.051527±0.000070" ]; then
    echo "Passed."
else
    echo "ERROR."
fi

echo -n "Test 3. Edge length mean and SD of cube.ply: "
# run the test
a=$($m -i input/cube.ply -edgeLength)
# check result
if [ "$a" == "edgeLength: 1.138071±0.000242" ]; then
    echo "Passed."
else
    echo "ERROR."
fi
