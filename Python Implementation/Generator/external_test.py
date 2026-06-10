import generator

df = generator.intialization()

###
#get config parameters for jiggle positions
config = generator.read_config("./Config/config_sample.yml")
rng = generator.set_rng_seed(config["seed"])
df_new =generator.jiggle_positions(df, config["rho"], config["sig2"], rng)

df.to_csv("synthetic_population_jiggled.csv", index=False)
df_new.to_csv("synthetic_population_jiggled_twice.csv", index=False)
