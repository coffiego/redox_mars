# This file is used to calc time responses of atmosphere to some changes

# Parameter settings are described in photochemistry_auto_run.jl
folder_directory_run = "/Users/shungokoyama/Programming/result/stat/"
including_file = folder_directory_run * "photochemistry_auto_run.jl"
include(including_file)
@everywhere @time update!(n_current,0.)

# we can run with logarithmic time steps
#@everywhere timepts = logspace(log10(1),log10(1e7*3.14e7),1000)
#@everywhere timediff = timepts[2:end]-timepts[1:end-1]
#@everywhere append!(timediff,3.3e6*3.14e7*ones(Float64,300))#koyama changed #originally 300, this runs for 2 billion years

@everywhere timepts = logspace(log10(1),log10(3.14e9),700)
@everywhere timepts_added = logspace(log10(3.14e9),log10(1e7*3.14e7),1600)
@everywhere append!(timepts,timepts_added)
@everywhere timediff = timepts[2:end]-timepts[1:end-1]
@everywhere append!(timediff,3.3e6*3.14e7*ones(Float64,300))

#run_filename="./CO2_20mbar_Tsurf_220_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2.h5"
core_filename = "CO2_" * string(CO2_pressure) * "mbar_Tsurf_" * string(surface_temperature) *
                "_Oesc_" *string(Oxygen_escape_rate_before)*"to"* string(Oxygen_escape_rate) *
                "_depv_" *string(depv) * "_nofixedCO2_noHOCO_Tmod.h5"
run_filename = folder_directory_run * core_filename
hfile = folder_directory_run*"Hesc_depfluxes_" * core_filename

@everywhere function runprofile(n_current, dtlist, filename)

    #このn_internalをこの関数の中で更新していく。
    n_internal = deepcopy(n_current)

    elapsed_time = 0.0

    # create a matrix to contain the data in n_internal, which is a Dict
    #それぞれのspeciesに対する高度プロファイル
    n_internal_mat = Array{Float64}(length(alt)-2,length(collect(keys(n_internal))));
    for ispecies in 1:length(collect(keys(n_internal)))
        for ialt in 1:length(alt)-2
            n_internal_mat[ialt,ispecies] = n_internal[collect(keys(n_internal))[ispecies]][ialt]
        end
    end

    #write initial input of the density profile
    h5write(filename,"n_current/init",n_internal_mat)
    h5write(filename,"n_current/alt",alt)
    h5write(filename,"n_current/species",map(string,collect(keys(n_internal))))
    h5write(filename,"n_current/timelist",cumsum(dtlist)) #cumulative sum along this vector

    # TIME LOOP - simulate over time
    # This section was changed so that it only prints statements as it starts
    # and finishes new tasks, rather than printing a statement at each dt.
    # the old lines for printing at each dt have been commented out and can be
    # restored if desired.
    thisi=0
    for dt in dtlist
        #println(filename*": iteration = "* string(thisi+=1)*" "*Libc.strftime(time()))
        thisi += 1
        #println("dt = "* string(dt::Float64))
        elapsed_time+=dt
        #println("elapsed_time = "*string(elapsed_time))
        # where the action happens - update n_internal for new timesteps
        update!(n_internal,dt)

        # save the concentrations to history
        # write n_current into n_current_mat
        n_internal_mat = Array{Float64}(length(alt)-2,length(collect(keys(n_internal))));
        for ispecies in 1:length(collect(keys(n_internal)))
            for ialt in 1:length(alt)-2
                n_internal_mat[ialt,ispecies] = n_internal[collect(keys(n_internal))[ispecies]][ialt]
            end
        end

        # write n_internal_mat to file
        h5write(filename, string("n_current/iter_",thisi), n_internal_mat)
    end
    return n_internal
end


@everywhere function read_ncurrent_from_file(readfile,tag)
    thisalt = h5read(readfile,"n_current/alt")
    if thisalt != alt
        throw("altitudes in file do not match altitudes in memory!")
    end
    n_current_tag_list = map(Symbol,h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,tag);
    n_current = Dict{Symbol,Array{Float64,1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]]=reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
end

@everywhere function get_H_fluxes(readfile)
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)

    Hfluxes=fill(0.,timelength)
    n_current=read_ncurrent_from_file(readfile,string("n_current/init"))
    Hfluxes[1]=(n_current[:H][end]*speciesbcs(:H)[2,2]
                  +2*n_current[:H2][end]*speciesbcs(:H2)[2,2])

    for i in 1:(timelength-1)
        n_current=read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        Hfluxes[i+1]=(n_current[:H][end]*speciesbcs(:H)[2,2]
                      +2*n_current[:H2][end]*speciesbcs(:H2)[2,2])
    end
    Hfluxes
end

#to plot deposition fluxes of O3,H2O2,HO2
@everywhere function get_all_outfluxes(readfile)
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)

    #create dominant deposition specise flux arrays
    Hfluxes=fill(0.,timelength)
    H2O2_dep=fill(0.,timelength)
    O3_dep=fill(0.,timelength)
    HO2_dep=fill(0.,timelength)
    H_dep=fill(0., timelength)
    totalH_outflux = fill(0.,timelength)
    totalO_outflux = fill(0.,timelength)

    n_current=read_ncurrent_from_file(readfile,string("n_current/init"))
    Hfluxes[1]=(n_current[:H][end]*speciesbcs(:H)[2,2]
                  +2*n_current[:H2][end]*speciesbcs(:H2)[2,2])
    H2O2_dep[1] = n_current[:H2O2][1]*speciesbclist[:H2O2][1,2]
    O3_dep[1] = n_current[:O3][1]*speciesbclist[:O3][1,2]
    HO2_dep[1] = n_current[:HO2][1]*speciesbclist[:HO2][1,2]
    H_dep[1] = n_current[:H][1]*speciesbclist[:H][1,2]
    totalH_outflux[1] = Hfluxes[1]+2*H2O2_dep[1]+HO2_dep[1]+H_dep[1] + n_current[:OH][1]*speciesbclist[:OH][1,2]
    totalO_outflux[1] = 2*H2O2_dep[1] + 3*O3_dep[1] + 2*HO2_dep[1] + speciesbclist[:O][2,2] + n_current[:OH][1]*speciesbclist[:OH][1,2] +n_current[:O][1]*speciesbclist[:O][1,2] +
                        n_current[:O1D][1]*speciesbclist[:O1D][1,2]

    for i in 1:(timelength-1)
        n_current=read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        Hfluxes[i+1]=(n_current[:H][end]*speciesbcs(:H)[2,2]
                      +2*n_current[:H2][end]*speciesbcs(:H2)[2,2])
        H2O2_dep[1+i] = n_current[:H2O2][1]*speciesbclist[:H2O2][1,2]
        O3_dep[1+i] = n_current[:O3][1]*speciesbclist[:O3][1,2]
        HO2_dep[1+i] = n_current[:HO2][1]*speciesbclist[:HO2][1,2]
        H_dep[1+i] = n_current[:H][1]*speciesbclist[:H][1,2]
        totalH_outflux[i+1] = Hfluxes[i+1]+2*H2O2_dep[i+1]+HO2_dep[i+1] + H_dep[i+1]+n_current[:OH][1]*speciesbclist[:OH][1,2]
        totalO_outflux[i+1] = 2*H2O2_dep[i+1] + 3*O3_dep[i+1] + 2*HO2_dep[i+1] + speciesbclist[:O][2,2] + n_current[:OH][1]*speciesbclist[:OH][1,2] +n_current[:O][1]*speciesbclist[:O][1,2] +
                            n_current[:O1D][1]*speciesbclist[:O1D][1,2]
    end
    (Hfluxes, H2O2_dep, O3_dep, HO2_dep, H_dep, totalH_outflux, totalO_outflux)
end

#ファイルを読み込んで, reactionratesとfluxes関数を使って入れて行くだけ
@everywhere function get_rates_and_fluxes(readfile)
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)
    reactionrateshist = fill(convert(Float64,NaN),timelength,length(intaltgrid),length(reactionnet))
    fluxhist = fill(convert(Float64,NaN),timelength,length(intaltgrid),length(specieslist))
    n_current = read_ncurrent_from_file(readfile,string("n_current/init"))
    reactionrateshist[1,:,:] = reactionrates(n_current)
    fluxhist[1,:,:] = fluxes(n_current,dz)
    for i in 1:(timelength-1)
        n_current = read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        reactionrateshist[i+1,:,:] = reactionrates(n_current)
        fluxhist[i+1,:,:] = fluxes(n_current,dz)
    end
    (reactionrateshist,fluxhist)
end

@everywhere function get_all_rates_and_fluxes(readfile)
    (reactionrateshist,fluxhist)=get_rates_and_fluxes(readfile)
    h5write(readfile,"fluxes/flux_history",fluxhist)
    h5write(readfile,"rates/reaction_rates_history",reactionrateshist)
    return
end

println("start to run!")
result = runprofile(n_current, timediff, run_filename)
println("finished running then start to make files for plot")
get_all_rates_and_fluxes(run_filename)


#This gets deposition get_H_fluxes
(Hfluxes,H2O2dep,O3dep,HO2dep,Hdep,totalH_outflux,totalO_outflux)=get_all_outfluxes(run_filename)
HOratio = totalH_outflux./totalO_outflux


# write out the H fluxes and deposition fluxes =======================================================
println("Writing escape and depostion fluxes file")
h5open(hfile, isfile(hfile) ? "r+" : "w") do file
   write(file,"fluxes/fluxvals",Hfluxes)
   write(file,"fluxes/times",h5read(run_filename,"n_current/timelist"))
   write(file, "depositions/H2O2",H2O2dep)
   write(file,"depositions/O3",O3dep)
   write(file,"depositions/HO2",HO2dep)
   write(file, "depositions/H", Hdep)
   write(file,"outfluxes/H",totalH_outflux)
   write(file,"outfluxes/O",totalO_outflux)
   write(file,"HOratio",HOratio)
   #write(file,"waterprofs/ppm",writewaterprof)
   #write(file,"waterprofs/alt",alt[2:end-1])
end

timelist = cumsum(timediff)
regtimeindex_rev = findfirst(x ->((x<1.99)|(2.01<x)), reverse(HOratio)) -1
regtimeindex = length(HOratio)-regtimeindex_rev + 1
println(run_filename)
if regtimeindex < 2
    println("not regulated")
elseif regtimeindex < length(timelist)
    println("regulation timescale: ", timelist[regtimeindex-1])
    h5open(hfile, isfile(hfile) ? "r+" : "w") do file
       write(file,"regulationtimescale",timelist[regtimeindex-1])
    end
else
    println("not regulated, need to check the result carefully")
end
#=
reg_header = ["CO2pressure","Ts","H2Oppm","Oesc","depv","regulationtimescale"]
reg_header = reshape(reg_header,1,length(reg_header))
open("summary_reg.csv","w") do io
    writedlm(io,reg_header,',')
end
=#
key_elements=[CO2_pressure,surface_temperature,H2Oppm,Oxygen_escape_rate_before,Oxygen_escape_rate,depv,timelist[regtimeindex-1]]
key_elements = reshape(key_elements,1,length(key_elements))
open(folder_directory_run*"/summary_reg_Tmod.csv","a") do io
    writedlm(io,key_elements,',')
end

println("done!")
