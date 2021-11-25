using Test, PlasmaProperties # This load both the test suite and our MyAwesomePackage
out = greet()
@test out == 123   