using Pkg
Pkg.add("BioStructures")
import BioStructures
using BioStructures
S
# Struct to track angle stats
global A = 0
global B = 0
global C = 0
global D = 0
global E = 0
global F = 0

# Function to print stats
function printStatus()
    println("A: $A || B: $B || C: $C || D: $D || E: $E || F: $F")
end

# Function to check if the angle is within a range
function checkrange(angle, upperbound, lowerbound)
    
    return (angle <= max(upperbound, lowerbound) && angle >= min(upperbound, lowerbound))
end


# Struct to hold pdb details

mutable struct pdb_details
    pdbCount
    pdbId
    resName_1
    resName_2
    resName_3
    resName_4
    resName_5
    resNumber_1
    resNumber_2
    resNumber_3
    resNumber_4
    resNumber_5
    distance_1
    distance_2
    distance_3
    distance_4
    psi_1
    psi_2
    psi_3
    phi_1
    phi_2
    phi_3
    turnType
    turnType_2
end






global pdb_count = 0

# Analysis function
function dataAnalysis(pdb_name::String)
    pdb_id = split(pdb_name, ".")[1]
    println("--------------", pdb_id, "--------------")
    global pdb_count += 1
    struc = read((joinpath(pdb_path, pdb_name)), PDBFormat) # reading pdb from files
    chains = BioStructures.collectchains(struc[1])

    for chain in chains
        for i in 1:(BioStructures.countresidues(chain))  # Avoid index out of bounds
            try
                res1 = chain[i]
                res2 = chain[i+1]
                res3 = chain[i+2]
                res4 = chain[i+3]
                res5 = chain[i+4]
                #res names for i to i+4
                r1 = resname(res1)
                r2 = resname(res2)
                r3 = resname(res3)
                r4 = resname(res4)
                r5 = resname(res5)
                #distances between O to N
                distance1 = BioStructures.distance(res1["O"], res4["N"])  # 2 < distance < 4.5
                distance2 = BioStructures.distance(res2["O"], res5["N"])
                #distances between C-alpha of i to C-alpha of i+3
                distance3 = BioStructures.distance(res1["CA"], res4["CA"]) # distance < 7.0
                distance4 = BioStructures.distance(res2["CA"], res5["CA"])

                #phi angles 
                phi1 = rad2deg(BioStructures.dihedralangle(res1["C"], res2["N"], res2["CA"], res2["C"]))
                phi2 = rad2deg(BioStructures.dihedralangle(res2["C"], res3["N"], res3["CA"], res3["C"]))
                phi3 = rad2deg(BioStructures.dihedralangle(res3["C"], res4["N"], res4["CA"], res4["C"]))

                #psi angles
                psi1 = rad2deg(BioStructures.psiangle(res2, res3))
                psi2 = rad2deg(BioStructures.psiangle(res3, res4))
                psi3 = rad2deg(BioStructures.psiangle(res4, res5))

                #default type if does not satisfy the psi and phi angle criteria
                type = "undefined"  #for single 
                type_2 = "undefined"  # for double

#single beta turn psi and phi angle criteria starts here
                if (checkrange(phi1, -30, -90) && checkrange(phi2, -60, -120) && checkrange(psi1, -1, -61) && checkrange(psi2, 45, -45))
                    type = "I"
                    global A += 1
                elseif (checkrange(phi1, 90, 30) && checkrange(phi2, 120, 60) && checkrange(psi1, 61, 1) && checkrange(psi2, 45, -45))
                    type = "I'"
                    global B += 1
                elseif (checkrange(phi1, -30, -90) && checkrange(phi2, 110, 50) && checkrange(psi1, 150, 90) && checkrange(psi2, 45, -45))
                    type = "II"
                    global C += 1
                elseif (checkrange(phi1, 90, 30) && checkrange(phi2, -50, -110) && checkrange(psi1, -90, -150) && checkrange(psi2, 45, -45))
                    type = "II'"
                    global D += 1
                elseif (checkrange(phi1, -30, -90) && checkrange(phi2, -30, -90) && checkrange(psi1, -1, -61) && checkrange(psi2, -1, -61))
                    type = "III"
                    global E += 1
                elseif (checkrange(phi1, 90, 30) && checkrange(phi2, 90, 30) && checkrange(psi1, 61, 1) && checkrange(psi2, 61, 1))
                    type = "III'"
                    global F += 1
                end
# remove this condition below for single beta turn
#double beta turn psi and phi angle criteria starts here
                if (checkrange(phi2, -30, -90) && checkrange(phi3, -60, -120) && checkrange(psi2, -1, -61) && checkrange(psi3, 45, -45))
                    type_2 = "I"
                    global A += 1
                elseif (checkrange(phi2, 90, 30) && checkrange(phi3, 120, 60) && checkrange(psi2, 61, 1) && checkrange(psi3, 45, -45))
                    type_2 = "I'"
                    global B += 1
                elseif (checkrange(phi2, -30, -90) && checkrange(phi3, 110, 50) && checkrange(psi2, 150, 90) && checkrange(psi3, 45, -45))
                    type_2 = "II"
                    global C += 1
                elseif (checkrange(phi2, 90, 30) && checkrange(phi3, -50, -110) && checkrange(psi2, -90, -150) && checkrange(psi3, 45, -45))
                    type = "II'"
                    global D += 1
                elseif (checkrange(phi2, -30, -90) && checkrange(phi3, -30, -90) && checkrange(psi2, -1, -61) && checkrange(psi3, -1, -61))
                    type_2 = "III"
                    global E += 1
                elseif (checkrange(phi2, 90, 30) && checkrange(phi3, 90, 30) && checkrange(psi2, 61, 1) && checkrange(psi3, 61, 1))
                    type_2 = "III'"
                    global F += 1
                end
#double beta turn condition ends

#distance criteria starts here
                if (checkrange(distance1, 4.5, 2.0) && checkrange(distance2, 4.5, 2.0) && distance3<7.0 && distance4<7.0 ) #for double beta turn
                # if  (checkrange(distance1, 4.5, 2.0 ) &&  distance3<7.0)    # for single beta turn
                try
                        p = pdb_details(pdb_count, pdb_id, r1, r2, r3, r4, r5, i, (i + 1), (i + 2), (i + 3), (i + 4), (round(distance1, digits=2)), (round(distance2, digits=2)), (round(distance3, digits=2)), (round(distance4, digits=2)), (Int(round(psi1))), (Int(round(psi2))), (Int(round(psi3))), (Int(round(phi1))), (Int(round(phi2))), (Int(round(phi3))), type ) #type_2 
                        printFileFunc(p)
                    catch
                    end
                end
            catch
            end
        end
    end
end




# File handling
global file_path = ""
global all_file = ""
global P_file = ""
global PG_file = ""
global PGA_file = ""


function openFiles()
    global file_path = raw"/Users/surajitmetya/Library/CloudStorage/OneDrive-Personal/Desktop/SouravMandal/" # paht to save data
    global all_file = open(joinpath(file_path, "all_single.csv"), "a")
    global undefined_file = open(joinpath(file_path, "undefine_single_pdb.csv"), "a")
    heading = "Count,ID,i,Res(i),i+1,Res(i+1),i+2,Res(i+2),i+3,Res(i+3),i+4,Res(i+4),D1,D2,D3,D4,Phi1,Phi2,Phi3,Psi1,Psi2,Psi3,Type1" #header for csv files
    println(all_file, heading)
    println(undefined_file, heading)
end

# Buffer to hold data
const buffer_size = 1000
global buffer_all = String[]
global buffer_undefined = String[]

# print in buffer process
function printFileFunc(p::pdb_details)
    data_to_print = "$(p.pdbCount),$(p.pdbId),$(p.resNumber_1),$(p.resName_1),$(p.resNumber_2),$(p.resName_2),$(p.resNumber_3),$(p.resName_3),$(p.resNumber_4),$(p.resName_4),$(p.resNumber_5),$(p.resName_5),$(p.distance_1),$(p.distance_2),$(p.distance_3),$(p.distance_4),$(p.phi_1),$(p.phi_2),$(p.phi_3),$(p.psi_1),$(p.psi_2),$(p.psi_3),$(p.turnType)," #$(p.turnType_2)
    
    if p.turnType != "undefined"
        push!(buffer_all, data_to_print)
        
    else
        push!(buffer_undefined, data_to_print)
    end
    
    # Flush buffers when  specified size is ready
    if length(buffer_all) >= buffer_size
        flushBuffers()
    end
end

# buffer to file --> then remove memory from buffer
function flushBuffers()
    write(all_file, join(buffer_all, "\n") * "\n")
    write(undefined_file, join(buffer_undefined, "\n") * "\n")

# Clear buffers
    empty!(buffer_all)
    empty!(buffer_undefined)
end

# Finalize writing and close all file handles
function finalizeWriting()
    flushBuffers()  # Flush any remaining data
    close(all_file)
    close(undefined_file)
end
pdb_path = raw"/Users/surajitmetya/Library/CloudStorage/OneDrive-Personal/Desktop/SouravMandal/PDB/All" # pdb files path
# Main function
@time function processPDBs()
    pdb_files =  [file for file in readdir(pdb_path) if endswith(file, ".pdb")] #list of pdb files
    openFiles()
    for pdb_file in pdb_files
        
        dataAnalysis(pdb_file)
    end
    finalizeWriting()
end

# Call the main function
processPDBs()
