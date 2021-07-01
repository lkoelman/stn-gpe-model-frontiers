#=
Spike train analysis functions

@author     Lucas Koelman
=#

using Distributions

"""
Function burst_metrics = legendy_new3(ISI,fac,sampling_rate,min_length_of_burst,local_length,surprise_cutoff)
carries out the burst detection algorithm by Legendy and Salcman (J. Neurophysiology 53:926) on 
the input ISI data stream.  The sampling rate defaults to 1000 Hz, but can be altered by giving
the sampling_rate as an input parameter.  As the initial criterion for burst detection, spikes 
have to occur at a frequency which is higher than the baseline frequency by a factor fac.  The
user can modify the minimal length of a prospective burst (min_length_of_burst), and the surprise
cutoff value.  In the original paper, this value was set to 10 to detect significant spikes.
The burst discharge is compared to a segment of data immediately preceding the burst.  A local_length
value of 0 indicates that the entire data stream should be used to yield the firing rate (can be used
in case of homogeneous firing of units).  Other values indicate a local time interval that should be
used (for instance, a value of 1 would mean that each spike/burst is compared to 1 second of data
prior to the index spike).
The function makes use of a subroutine called surprise.  The function returns 100 for
a very high surprise value and 0 for a very low one.  This means that the output of this rou-
tine is not accurate for such very high or very low surprise values (although bursts are correctly
detected).  
The function produces a structure (burst_metrics) which contains fields describing the bursts, including the 
onset and lengths of bursts (described as indices of the input ISI stream) the highest rate within 
the burst, the average discharge rate in the burst, the pre-burst baseline rate, as well as the 
surprise values for individual bursts.  In field 1 of the structure, summary parameters are defined, 
including, num_bursts (the total number of bursts detected), mean_spikes_per_burst, total_spikes_in_bursts, 
mean_intra_burst_frequency, proportion_time_in_bursts, and proportion_spikes_in_bursts.
Defaults: surprise_cutoff = 10, local_length = 0, min_length_of_burst = 2, sampling_rate = 1000, fac = 1.5
Usage: burst_metricss = legendy_new3(A,2,10000,4,1000,10)
This call would detect bursts in the input ISI stream A, sampled at 10000 Hz ISIs within such bursts
would have to fire at at least twice the local firing rate, and would have to be at least 4 spikes 
in length.  The comparison firing rates would be computed using 1000 spikes preceding the prospective 
bursts.  burst_metricss would have to have a surprise value of at least 10.

Written 9/21-23/2001, 7/24/2003, 11/27/2003, 7/31/2005, 1/4/2007 by Thomas Wichmann.
"""
function burst_metrics_surprise(
    ISI,
    fac=1.5,
    sampling_rate=1000,
    min_length_of_burst=2,
    local_length=0,
    surprise_cutoff=10)

burst_metrics = Dict(
        "begin" => [],
        "num_spikes"=> [],
        "surprise" => [],
        "rate" => [],
        "max_rate"=> [],
        "baseline_rate" => [],
        # Summary metrics
        "num_bursts" => 0,
        "mean_spikes_per_burst" => 0,
        "median_spikes_per_burst" => 0,
        "total_spikes_in_bursts" => 0,
        "mean_intra_burst_frequency" => 0.0,
        "median_intra_burst_frequency" => 0.0,
        "proportion_time_in_bursts" => 0.0,
        "proportion_spikes_in_bursts" => 0.0
    )

burst_num = 0
CA = cumsum(ISI)

if local_length == 0
    mean_FR = length(ISI) / (sum(ISI)/sampling_rate)
    fr_thr = sampling_rate / (fac * mean_FR)   # calculation of frequency threshold
    beg_idx = 0
else
    # finds last index within the 'local length' - incremented by one, this will result in the first index that can be evaluate.
    beg_idx = findlast(CA .< local_length*sampling_rate)
end
n = beg_idx

# ***** Main loop ****

while n < length(ISI) - min_length_of_burst

    n = n+1 # n is a running parameter that points to ISIs
    
    if local_length > 0
        # find the ISI data segment I that is fully contained within the local_length
        I = ISI[findfirst(CA .> CA(n)-local_length*sampling_rate):n-1]
        mean_FR = length(I) / (sum(I) / sampling_rate)
        fr_thr = sampling_rate / (fac*mean_FR)  # calculation of frequency threshold
    end
    
    
    # ****** 1. selection step - find runs of short ISIs *******************
    if (ISI[n] < fr_thr) # finding areas of the spike train which fulfill the length_of_burst criterion

        q = 0           # running parameter that points to the number of spikes to be added
        inc_flag = false
        while (n+q <= length(ISI)) && (ISI[n+q] < fr_thr)
            q = q + 1
            inc_flag = true
        end
        if inc_flag
            q = q - 1                                                                # reverse the last addition of q that led to the end of the while routine
        end
        
        # at this point, the provisional burst starts at n and ends at n+q
        # it has q+1 spikes in it
    
        
    # ******* 2. selection step - adjust length of burst to maximize surprise value ***************
        if q+1 >= min_length_of_burst
            m = min_length_of_burst                                            # running parameter of the number of spikes to be added to n
            inc_flag = false
            while ((n+m <= length(ISI)) &&
                    (ISI[n+m] < fr_thr) &&
                    (surprise(mean_FR, ISI[n:n+m], sampling_rate) >= 
                     surprise(mean_FR, ISI[n:n+m-1], sampling_rate)))   # 2. burst criterion - surprise values are increasing when spikes are added 
                m = m+1
                inc_flag = true
            end
            if inc_flag
                m = m-1                            # reverse the last addition steps that resulted in the termination of the while loop
            end
           
            # at this point, the beginning of the burst is provisionally settled to be n, the end n+m
            # the burst has m+1 spikes in it.
            
    # ******* 3. selection step - test whether adding up to 10 more ISIs will enhance surprise value **********************
            if n+m+10 <= length(ISI) # mmax is set to 10 unless one couldn't add 10 to n before reaching the end of FR
                mmax = 10
            else
                mmax = length(ISI)-(n+m)
            end
        
            ind_long_ISI = findfirst(ISI[n+m+1:n+m+mmax] .> fr_thr)        # looking for 'slow spikes' within the next 10 spikes after the current burst end
            if ~isempty(ind_long_ISI)                                           # pmax is set to be the index of the slow spike
                pmax = ind_long_ISI-1
            else
                pmax = mmax
            end
            
            S = zeros(pmax+1)                                                # formation of an array that will contain surprise values.  The first one is that of the burst defined by the ISIs between n and n+m.  Additional entries into S are surprise values that would result from adding up to pmax additional spikes (S2 will be the surprise values for ISI[n:n+m+1], S3 the one for ISI[n:n+m+2] etc.)
            S[1] = surprise(mean_FR, ISI[n:n+m], sampling_rate)                                        
            for p = 1:pmax                                  # forms array of surprise values for this burst, starting from the end of the burst to pmax values later
                S[p+1] = surprise(mean_FR, ISI[n:n+m+p], sampling_rate)
            end
            ind_max_S = argmax(S)            

            if n+m < length(ISI)
                m = m + ind_max_S - 1                          # this will set the m-value to the previous m, if the first entry into the S array is the highest (i.e., if adding additional ISIs didn't increase the surprise value), or, it will correct m to the highest value
            else
                m = length(ISI)-n
            end
        
            # at this point, the end of the index of the end of the burst
            # is settled to be n+m
            
            
        # ******** 4. selection step - test whether adjusting the front end of the burst enhances the surprise value ******************
            if n > 1
                o = 1 
                inc_flag = false
                while((m-o > min_length_of_burst) &&
                      (surprise(mean_FR, ISI[n+o:n+m], sampling_rate) >= 
                       surprise(mean_FR, ISI[n+o-1:n+m], sampling_rate)))
                    o = o+1
                    inc_flag = true
                end
                
                if inc_flag
                    o = o - 1          # reducing o by one to correct for the addition that resulted in the end of the while loop
                    n = n+o            # adjust the beginning of the burst
                    m = m-o            # adjust the length of the burst
                end
            end
        
            # at this point, the beginning of the burst is settled to be n, and the length is m+1 ***
            
            if (m+1 >= min_length_of_burst) && (surprise(mean_FR,ISI[n:n+m],sampling_rate) > surprise_cutoff)
                
                burst_num = burst_num + 1
                push!(burst_metrics["begin"], n)
                push!(burst_metrics["num_spikes"], m+1)
                push!(burst_metrics["surprise"], surprise(mean_FR,ISI[n:n+m], sampling_rate))
                push!(burst_metrics["rate"], length(ISI[n:n+m]) / (sum(ISI[n:n+m])/sampling_rate))
                push!(burst_metrics["max_rate"], sampling_rate/minimum(ISI[n:n+m]))
                push!(burst_metrics["baseline_rate"], mean_FR)
            end        
            
            n = n+m+1  # adjust ISI pointer to the ISI following the burst
            
        end
    end
end


# ****** Store burst parameters in the output array
if !isempty(burst_metrics["begin"])

    num_bursts = length(burst_metrics["begin"])
    burst_metrics["num_bursts"] = num_bursts

    C = zeros(num_bursts)
    for n = 1:num_bursts
        C[n] = burst_metrics["num_spikes"][n]
    end
    burst_metrics["mean_spikes_per_burst"] = mean(C)
    burst_metrics["median_spikes_per_burst"] = median(C)
    burst_metrics["total_spikes_in_bursts"] = sum(C)
    
    for n = 1:num_bursts
        C[n] = burst_metrics["rate"][n]
    end
    burst_metrics["mean_intra_burst_frequency"] = mean(C)
    burst_metrics["median_intra_burst_frequency"] = median(C)
    burst_metrics["proportion_time_in_bursts"] = (burst_metrics["total_spikes_in_bursts"] / burst_metrics["mean_intra_burst_frequency"]) / (sum(ISI[beg_idx+1:length(ISI)]) / sampling_rate)
    
    burst_metrics["proportion_spikes_in_bursts"] = burst_metrics["total_spikes_in_bursts"] / length(ISI[beg_idx+1:length(ISI)])
else 
    burst_metrics["num_bursts"] = 0
    burst_metrics["mean_spikes_per_burst"] = 0
    burst_metrics["median_spikes_per_burst"] = 0
    burst_metrics["total_spikes_in_bursts"] = 0
    burst_metrics["mean_intra_burst_frequency"] = 0
    burst_metrics["median_intra_burst_frequency"] = 0
    burst_metrics["proportion_time_in_bursts"] = 0
    burst_metrics["proportion_spikes_in_bursts"] = 0
end

return burst_metrics

end # end function

#**********************************************************************

"""
Calculates surprise index.

# Arguments

r = comparison firing rate (spikes per second)
data = ISI data to be included in the burst
sampling_rate = sampling rate of ISI measurements
"""
function surprise(r, data, sampling_rate; with_deceleration=false)


    T = sum(data) / sampling_rate
    num_spikes = length(data)

    poiss = Poisson(r*T)
    p = cdf(poiss, num_spikes)

    if p == 0           # for very small p, a value of 10exp-100 is assumed
        burst = 0
        deceleration = 100
    elseif p == 1      # for very high p, a value of 1-10exp-100 is assumed    
        burst = 100    
        deceleration = 0
    else
        burst = -log(1-p)
        deceleration = -log10(p)
    end

    if with_deceleration
        return (burst, deceleration)
    else
        return burst
    end

end # end function


"""
Test for burst_metrics_surprise()
"""
function test_burst_metrics()

    times = collect(1:1000)
    jitter = rand(Float64, length(times))
    spikes = sort(times + jitter * 20.0)
    isi_values = diff(spikes)

    metrics = burst_metrics_surprise(isi_values)

end