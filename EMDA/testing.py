from emda import EMDA

emda = EMDA("../example/parameters.prmtop", "../example/trajectory.nc")
emda.select("first_resids", [1, 2, 3], sel_type="res_num")
emda.select("second_resids", [4, 5, 6], sel_type="res_num")

emda.add_distance(
    name="dist_first_second", sel1="first_resids", sel2="second_resids", type="min"
)
emda.run()

print(emda.measures["dist_first_second"].result)
