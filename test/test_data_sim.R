data_settings_in <- list(lod=45,
                         n=300,
                         p_group=c(0.3,0.6,0.2),
                         t_obs=(-6):7,
                         sensitivity=0.99)

param_settings_in <- list(p_model=c(0.3,0.5,0.2),
                          wp_mean=c(5,5,5),
                          wp_sd=c(2,2,2),
                          tp_sd=c(1,1,1),
                          dp_mean=c(30,30,30),
                          dp_sd=c(5,5,5),
                          wr_mean=c(7,7,7),
                          wr_sd=c(2,2,2),
                          sigma=5,
                          wp_min=1,
                          wp_max=8,
                          dp_min=20,
                          tp_min=-1.5,
                          tp_max=1.5,
                          wr_min=1,
                          wr_max=10)

sim_out <- simulate_viral_load_data(data_settings_in,
                                    param_settings_in,
                                    sim_seed=12)

