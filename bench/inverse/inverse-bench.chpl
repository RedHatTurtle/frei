config const n : int = 3;
config const testTime : real = 60;

proc main()
{
    use Time;
    use Random;
    use LinearAlgebra;
    use Inverse;

    var matrix : [0..<n, 0..<n] real;

    var matrixInv : [0..<n, 0..<n] real;

    // Generate random test matrix
    fillRandom(matrix, 17);

    writeln("Test duration: ", testTime);

    var watch : stopwatch;
    var counter : int = 0;
    watch.start();
    while watch.elapsed() < testTime
    {
        matrixInv = inv(M=matrix);
        counter += 1;
    }

    writef("MatrixInv     reps: %i,%03i\n", counter/1000, counter%1000);

    counter = 0;
    watch.restart();
    while watch.elapsed() < testTime
    {
        matrixInv = inv_dot(M=matrix);
        counter += 1;
    }

    writef("MatrixInv_dot reps: %i,%03i\n", counter/1000, counter%1000);
    writeln();
}
