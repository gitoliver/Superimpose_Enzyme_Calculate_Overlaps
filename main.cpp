#include <fstream>
#include <iostream>
#include <cstdlib> // for exit()
// Includes for directory reading
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>     /* getenv */
#include <thread>
#include <future>
#include <sys/wait.h>
#include "io.h"
#include "residue_name_conversion_map.h"
#include "/home/oliver/Programs/gems/gmml/includes/gmml.hpp"


typedef std::vector<std::string> StringVector;
typedef std::vector<Assembly*> AssemblyVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<MolecularModeling::Atom*> AtomVector;


constexpr auto PI = 3.14159265358979323846;
constexpr auto CSA = 36.31681103; // Carbon Surface Area for normalizing results to the area of a carbon atom. Can ask the question: How many C atom equivalents are inside the protein?

Residue* find_residue_by_number(std::string residue_number, Assembly *assembly);
Residue* seek_alycone_residue_connected_to_residue_number(std::string residue_number, Assembly *assembly);
Residue* recursively_seek_alycone_residue(Atom *start_atom, bool *found);
double CalculateAtomicOverlapsForAssemblyVector(Assembly *assemblyA, AssemblyVector assemblies);
bool CheckIfProtein(std::string resname);
double SuperimposeSnapshotsToEnzymeCalculatePercentTimeAccessible(Assembly *enzyme, Assembly *enzyme_ligand, AssemblyVector substrate_snapshots, std::string substrate_residue_number);
//double SuperimposeEnzymeCalculateOverlaps(Assembly *enzyme, AtomVector target_atoms, AtomVector substrate, AtomVector substrate_superimposition_atoms);
AtomVector GetAtomSelection(Assembly *assembly, StringVector residue_numbers, StringVector atom_names);

int main(int argc, char *argv[])
{
    std::string working_Directory = Find_Program_Working_Directory();
    std::string installation_Directory = Find_Program_Installation_Directory();


    //************************************************//
    // Reading input file                             //
    //************************************************//

    std::string parameter_directory, enzyme_PDB, enzyme_ligand_PDB, substrate_directory, substrate_residue_numbers, residue_conversion_map;
    std::ifstream inf (working_Directory + "/inputs/" + "input.txt");
    if (!inf)
    {
        std::cerr << "Uh oh, input file could not be opened for reading!" << std::endl;
        std::exit(1);
    }
    while (inf) // While there's still stuff left to read
    {
        std::string strInput;
        getline(inf, strInput);
        if(strInput == "Parameters:")
            getline(inf, parameter_directory); // Reads a single line after a line with the keyword "Parameters:"
        if(strInput == "Enzyme:")
            getline(inf, enzyme_PDB);
        if(strInput == "Enzyme ligand:")
            getline(inf, enzyme_ligand_PDB);
        if(strInput == "Substrate directory:")
            getline(inf, substrate_directory);
        if(strInput == "Substrate residue numbers:")
            getline(inf, substrate_residue_numbers);
        if(strInput == "Residue Conversion Map:")
            getline(inf, residue_conversion_map);
    }

    // Map to convert from NLN residue numbers in input file to HXB2 numbering via user provided map.
    std::string map_input_file = (working_Directory + "/inputs/" + residue_conversion_map);
    Residue_Name_Conversion_Map conversion_map(map_input_file);
    conversion_map.PrintMap();

    if (fileExists(parameter_directory + "/amino12.lib"))
        std::cout << "Using user provided parameters" << std::endl;
    else
        std::cout << "Using default parameters from installation folder" << std::endl;
        parameter_directory = installation_Directory + "/CurrentParams";

    //************************************************//
    // Details for loading in a PDB file              //
    //************************************************//

    std::vector<std::string> amino_libs, glycam_libs, other_libs, prep;
    amino_libs.push_back(parameter_directory + "/amino12.lib");
    amino_libs.push_back(parameter_directory + "/aminoct12.lib");
    amino_libs.push_back(parameter_directory + "/aminont12.lib");

    glycam_libs.push_back(parameter_directory + "/GLYCAM_amino_06j_12SB.lib");
    glycam_libs.push_back(parameter_directory + "/GLYCAM_aminoct_06j_12SB.lib");
    glycam_libs.push_back(parameter_directory + "/GLYCAM_aminont_06j_12SB.lib");

    other_libs.push_back(parameter_directory + "/nucleic12.lib");
    other_libs.push_back(parameter_directory + "/nucleic12.lib");
    other_libs.push_back(parameter_directory + "/solvents.lib");

    prep.push_back(parameter_directory + "/GLYCAM_06j-1.prep");

    //std::string pdb_file_path = working_Directory + "/outputs/" + "test.pdb";
    std::string parameter_file_path = parameter_directory + "/GLYCAM_06j.dat";
   // std::string ion_parameter_file_path = parameter_directory + "/atomic_ions.lib";

    //************************************************//
    // Load PDB files                                 //
    //************************************************//

    Assembly enzyme;
    enzyme.BuildAssemblyFromPdbFile( (working_Directory + "/inputs/" + enzyme_PDB), amino_libs, glycam_libs, other_libs, prep, parameter_file_path );
    enzyme.BuildStructureByDistance();
    AtomVector neighbors =enzyme.GetAllAtomsOfAssembly().at(0)->GetNode()->GetNodeNeighbors();
    std::cout << "Number of neighbors is " << neighbors.size() << std::endl;

    //AtomVector test = enzyme.Select(":#4");
    //std::cout << "Selected " << test.at(0)->GetResidue()->GetId() << std::endl;

    Assembly enzyme_ligand;
    enzyme_ligand.BuildAssemblyFromPdbFile( (working_Directory + "/inputs/" + enzyme_ligand_PDB), amino_libs, glycam_libs, other_libs, prep, parameter_file_path );
    enzyme_ligand.BuildStructureByDistance();

    //************************************************//
    // Load files from directory                      //
    //************************************************//

    std::string directory = working_Directory + "/inputs/" + substrate_directory ;
    std::cout << "directory: " << directory << std::endl;
    std::string filepath;
    AssemblyVector substrate_snapshots;

    DIR *dp; // A directory stream
    struct dirent *dirp; // Contains file serial number and name (char d_name[])
    struct stat filestat; // Contains info about file, such as device ID, user ID, access time etc

    dp = opendir( directory.c_str() ); //.c_str adds a null character to the end.
    if (dp == NULL)
    {
        std::cout << "Error(" << errno << ") opening " << directory << std::endl;
        return errno;
    }
    while ((dirp = readdir ( dp )))
    {
        filepath = directory + "/" + dirp->d_name;
        // If the file is a directory (or is in some way invalid) we'll skip it
        if (stat( filepath.c_str(), &filestat )) continue; // Is it a valid file?
        if (S_ISDIR( filestat.st_mode ))         continue; // Is it a directory?
        StringVector temp;
        temp.push_back(filepath);
        substrate_snapshots.push_back(new Assembly(temp, gmml::InputFileType::PDB));
        //temp_assembly.BuildAssemblyFromPdbFile(filepath, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
        //temp_assembly.BuildStructureByDistance();
    }
    closedir( dp );


    //************************************************//
    // Superimposition                                //
    //************************************************//
    std::cout << "Superimposition and overlap calculation" << std::endl;


/*
    //AtomVector alsoMoving_atoms = enzyme.GetAllAtomsOfAssembly();
    AtomVector target_atoms;
    AtomVector enzyme_ligand_atoms = enzyme_ligand.GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = enzyme_ligand_atoms.begin(); it != enzyme_ligand_atoms.end(); ++it)
    {
        Atom *atom = *it;
        if ( (atom->GetName().compare("C1")==0) || (atom->GetName().compare("C3")==0) || (atom->GetName().compare("C5")==0) )
            target_atoms.push_back(atom);
    }
    */


    std::string converted_resnum;
    substrate_snapshots.at(0)->BuildStructureByDistance();

    StringVector split_substrate_residue_numbers = split(substrate_residue_numbers, ',');
    for(StringVector::iterator it = split_substrate_residue_numbers.begin(); it != split_substrate_residue_numbers.end(); ++it)
    {
       // std::cout << "id is " << *it << std::endl;

        double percent_accessible = SuperimposeSnapshotsToEnzymeCalculatePercentTimeAccessible(&enzyme, &enzyme_ligand, substrate_snapshots, *it);

        Residue* algycone = seek_alycone_residue_connected_to_residue_number(*it, substrate_snapshots.at(0));
        converted_resnum = conversion_map.ConvertViaMap(algycone->GetNumber());
        std::cout << "Time accessible for " << converted_resnum << ": " <<  percent_accessible << std::endl;
    }
/*
   // int i=0; //counter only for output pdb file names
    for(AssemblyVector::iterator it = substrate_snapshots.begin(); it != substrate_snapshots.end(); ++it)
    {     
        Assembly *assembly = *it;
        AtomVector alsoMoving_atoms = assembly->GetAllAtomsOfAssembly();
        ResidueVector residues = assembly->GetAllResiduesOfAssembly();
        AtomVector super_atoms;
        for(ResidueVector::iterator itt = residues.begin(); itt != residues.end(); ++itt)
        {
            Residue *residue = *itt;
            StringVector split_id = split(residue->GetId(), '_');
            if ( substrate_residue_number.compare(split_id.at(2))==0 )
            {
                AtomVector atoms = residue->GetAtoms();
                for(AtomVector::iterator ittt = atoms.begin(); ittt != atoms.end(); ++ittt)
                {
                    Atom *atom = *ittt;
                    if ( (atom->GetName().compare("C1")==0) || (atom->GetName().compare("C3")==0) || (atom->GetName().compare("C5")==0) )
                        super_atoms.push_back(atom);
                }
            }
        }

        gmml::Superimpose(super_atoms, target_atoms, alsoMoving_atoms);

        //Write out a pdb file:
       // std::stringstream ss;
       // ss << working_Directory + "/outputs/moved_" << i << ".pdb";
       // PdbFileSpace::PdbFile *outputPdbFile = assembly->BuildPdbFileStructureFromAssembly(-1,0);
       // outputPdbFile->Write(ss.str());
      //  ++i;

       // double overlaps = CalculateAtomicOverlaps(&enzyme, assembly);
       // std::cout << CalculateAtomicOverlaps(&enzyme, assembly) << std::endl;
       // std::ofstream myfile;

      //  std::stringstream ss1;
      //  ss1 << working_Directory + "/outputs/Overlaps.txt";
     //   std::cout << ss1.str() << std::endl;
     //   myfile.open (ss1.str(), std::ios::out | std::ios::app);
     //   std::string filename = assembly->GetName() + ": ";
     //   myfile << filename;
     //   myfile << overlaps << "\n";
     //   myfile.close();
    }

    CalculateAtomicOverlapsForAssemblyVector(&enzyme, substrate_snapshots);
*/
    return 0;
}


/*bool CheckIfProtein(std::string resname){
    std::vector<std::string> residue_list = {"ALA","ASP", "ASN", "ARG", "GLY", "GLU", "GLN", "PRO", "HIS", "CYS", "VAL", "LEU", "THR", "SER", "LYS", "MET", "TYR", "TRP", "PHE", "SEC", "ILE", "CYX", "HID", "HIE" };
    //std::string resname = this->GetName();
    if (std::any_of(residue_list.begin(), residue_list.end(),
                [resname](std::string residue_names) {
        if (residue_names.compare(resname)==0)
            return true;
        return false;
    }))
        return true;
    return false;
}*/

double SuperimposeSnapshotsToEnzymeCalculatePercentTimeAccessible(Assembly *enzyme, Assembly *enzyme_ligand, AssemblyVector substrate_snapshots, std::string substrate_residue_number)
{

    //Strings for selection of atoms for superimposition
    StringVector residue_selection_strings;
    residue_selection_strings.push_back(substrate_residue_number);
    StringVector atom_selection_strings = {"C1", "C3", "C5"};

    AtomVector enzyme_ligand_atoms = enzyme_ligand->GetAllAtomsOfAssembly();
    AtomVector target_atoms;
    for(StringVector::iterator it = atom_selection_strings.begin(); it != atom_selection_strings.end(); ++it)
    {
        for(AtomVector::iterator it1 = enzyme_ligand_atoms.begin(); it1 != enzyme_ligand_atoms.end(); ++it1)
        {
            Atom *atom = *it1;
            if (atom->GetName().compare(*it)==0)
            {
                target_atoms.push_back(atom);
            }
        }
    }

    int snapshot_accessible_counter = 0;

    // This next part is for running on multiple CPUs
    pid_t pids[substrate_snapshots.size()];
    int i;
   // unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
   // std::cout << "Number of cores: " << concurentThreadsSupported << std::endl; // Not using now, try optimize with this in future.
    int n = substrate_snapshots.size();

    std::cout << "Starting " << n << " children" << std::endl;
    /* Start children. */
    for (i = 0; i < n; ++i)
    {
        if ((pids[i] = fork()) < 0)
        {
            perror("fork");
            abort();
        }
        else if (pids[i] == 0)
        {
            //DoWorkInChild();
            AtomVector super_atoms = GetAtomSelection(substrate_snapshots.at(i), residue_selection_strings, atom_selection_strings );
            AtomVector substrate_all = substrate_snapshots.at(i)->GetAllAtomsOfAssembly();
            AtomVector substrate_protein = substrate_snapshots.at(i)->GetAllAtomsOfAssemblyWithinProteinResidues();
            //std::cout << "Substrate_protein.size " << substrate_protein.size() << std::endl;
           // substrate_snapshots.at(i)->BuildStructureByDistance(4);
            //substrate_snapshots.at(i)->BuildStructureByDistanceByOptimizedThread();
           // StringVector options = { "cutoff = 1.7",
          //  substrate_snapshots.at(i)->BuildStructure(gmml::DISTANCE);
           // std::cout << "X value is then " << substrate_protein.at(0)->GetCoordinates().at(0)->GetX() << std::endl;
            gmml::Superimpose(super_atoms, target_atoms, substrate_all);
          //  std::cout << "X value is now " << substrate_protein.at(0)->GetCoordinates().at(0)->GetX() << std::endl;
          //  std::cout << "Original is now " << substrate_snapshots.at(i)->GetAllAtomsOfAssembly().at(0)->GetCoordinates().at(0)->GetX() << std::endl;

            substrate_snapshots.at(i)->BuildStructureByDistance(4);
            std::stringstream ss;
            ss << "outputs/moved_" << i << "_" << substrate_residue_number << ".pdb";
            PdbFileSpace::PdbFile *outputPdbFile = substrate_snapshots.at(i)->BuildPdbFileStructureFromAssembly(-1,0);
            outputPdbFile->Write(ss.str());


            double current_overlap = enzyme->CalculateAtomicOverlaps(substrate_protein);
            std::cout << "enzyme overlap with substrate is " << current_overlap << std::endl;
            //double current_overlap = SuperimposeEnzymeCalculateOverlaps(enzyme, target_atoms, substrate, super_atoms);
            if (current_overlap < 2.5 )
            {
                ++snapshot_accessible_counter;
            }
            exit(0);
        }
    }

    /* Wait for children to exit. */
    int status;
    pid_t pid;
    while (n > 0)
    {
        pid = wait(&status);
        //  printf("Child with PID %ld exited with status 0x%x.\n", (long)pid, status);
        --n;  // TODO(pts): Remove pid from the pids array.
    }

    /*i = 0;
    for(AssemblyVector::iterator it = substrate_snapshots.begin(); it != substrate_snapshots.end(); ++it)
    {
        Assembly *substrate = *it;
        std::cout << "X value will now " << substrate->GetAllAtomsOfAssembly().at(0)->GetCoordinates().at(0)->GetX() << std::endl;
        //thread this:

       // futures_returns.push_back(std::async(&SuperimposeEnzymeCalculateOverlaps, enzyme, target_atoms, substrate, substrate_residue_number));


        //   double current_overlap = ret.get();
        //if (ret.get() < 2.5) { ++snapshot_accessible_counter; }
        // std::cout << "enzyme overlap with substrate is " << ret.get() << std::endl;
        //double current_overlap = SuperimposeEnzymeCalculateOverlaps(enzyme, target_atoms, substrate, substrate_residue_number);
        //  std::cout << "enzyme overlap with substrate is " << current_overlap << std::endl;

     //   if (current_overlap < 2.5 )
     //   {
     //       ++snapshot_accessible_counter;
    //    }
        //Write out a pdb file:
        substrate->BuildStructureByDistance(4);
        std::stringstream ss;
        ss << "outputs/moved_" << i << ".pdb";
        PdbFileSpace::PdbFile *outputPdbFile = substrate->BuildPdbFileStructureFromAssembly(-1,0);
        outputPdbFile->Write(ss.str());
        ++i;
    }*/

    if ( (snapshot_accessible_counter > 0) && (substrate_snapshots.size() > 0) )
        return ( (snapshot_accessible_counter/substrate_snapshots.size() ) * 100 );
    else
        return 0;
}

AtomVector GetAtomSelection(Assembly *assembly, StringVector residue_numbers, StringVector atom_names)
{
    ResidueVector residues = assembly->GetAllResiduesOfAssembly();
    ResidueVector selectedResidues;
    AtomVector selection;
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); ++it)
    {
        Residue *residue = *it;

        StringVector split_id = split(residue->GetId(), '_');
        //std::cout << "resid: " << split_id.at(2) << std::endl;
        for(StringVector::iterator it1 = residue_numbers.begin(); it1 != residue_numbers.end(); ++it1)
        {
            //std::string *resnum = *it1;
            if (it1->compare(split_id.at(2))==0)
            {
                selectedResidues.push_back(residue);
            }
        }
    }
    for(ResidueVector::iterator it = selectedResidues.begin(); it!= selectedResidues.end(); ++it)
    {
        Residue *residue = *it;
        AtomVector residueAtoms = residue->GetAtoms();
        for (StringVector::iterator it1 = atom_names.begin(); it1 != atom_names.end(); ++it1)
        {
            for (AtomVector::iterator it2 = residueAtoms.begin(); it2!= residueAtoms.end(); ++it2)
            {
                Atom *atom = *it2;
                if (it1->compare(atom->GetName())==0)
                {
                    //std::cout << "Adding " << residue->GetName() << ":" << atom->GetName() << " to selection" << std::endl;
                    selection.push_back(atom);
                }
            }
        }
    }
    return selection;
}

MolecularModeling::Residue* find_residue_by_number(std::string residue_number, Assembly *assembly)
{
    ResidueVector residues = assembly->GetAllResiduesOfAssembly();
    StringVector split_id;
    Residue *residue = new Residue();
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); ++it)
    {
        residue = *it;
        split_id = split(residue->GetId(), '_');
        if (residue_number.compare(split_id.at(2))==0)
        {
            return residue;
        }
    }
    std::cout << "Residue with number " << residue_number << " not found!" << std::endl;
    return residue;
}


MolecularModeling::Residue* seek_alycone_residue_connected_to_residue_number(std::string residue_number, Assembly *assembly)
{
    bool found = false;
    bool *pfound = &found;
    MolecularModeling::Residue *current_residue = find_residue_by_number(residue_number, assembly);
    Atom *start_atom = current_residue->GetAtoms().at(0); // Just get any atom to start from.
    Residue *aglycone = recursively_seek_alycone_residue(start_atom, pfound);
    //std::cout << "Out a bit with the aglycone: " << aglycone->GetId() << std::endl;
    return aglycone;
}

MolecularModeling::Residue* recursively_seek_alycone_residue(Atom *start_atom, bool *found)
{
    Residue *aglycone = new Residue();
    start_atom->SetDescription("Visited"); // Hmm. Not sure this is a good idea, but it sure is handy.
    //std::cout << "Checking neighbors of " << start_atom->GetId() << std::endl;
    AtomVector neighbors = start_atom->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
    {
        // If haven't found aglycone yet and haven't visited this neighbor already
        if ( (*found == false) && ((*it)->GetDescription().compare("Visited")!=0) )
        {
            Atom *neighbor = *it;
            if(neighbor->GetResidue()->CheckIfProtein())
            {
                aglycone = neighbor->GetResidue();
                *found = true;
            }
            else
            {
                aglycone = recursively_seek_alycone_residue(neighbor, found);
            }
        }
    }
    return aglycone;
}

//double SuperimposeEnzymeCalculateOverlaps(Assembly *enzyme, AtomVector target_atoms, AtomVector alsoMoving_atoms, AtomVector substrate_superimposition_atoms)
//{
    //std::cout << "Just launched superimposeThread" << std::endl;
    //AtomVector alsoMoving_atoms = substrate->GetAllAtomsOfAssembly();
    //ResidueVector residues = substrate->GetAllResiduesOfAssembly();
    //AtomVector super_atoms;
    /*
    for(ResidueVector::iterator itt = residues.begin(); itt != residues.end(); ++itt)
    {
        Residue *residue = *itt;
        StringVector split_id = split(residue->GetId(), '_');
        //std::cout << "resid: " << split_id.at(2) << std::endl;
        if ( substrate_residue_number.compare(split_id.at(2))==0 )
        {
            AtomVector atoms = residue->GetAtoms();
            for(AtomVector::iterator ittt = atoms.begin(); ittt != atoms.end(); ++ittt)
            {
                Atom *atom = *ittt;
                if ( (atom->GetName().compare("C1")==0) || (atom->GetName().compare("C3")==0) || (atom->GetName().compare("C5")==0) )
                {
                    super_atoms.push_back(atom);
                }
            }
        }
    }
    */

   // gmml::Superimpose(substrate_superimposition_atoms, target_atoms, alsoMoving_atoms);
    // Pull out closest X atoms from substrate here. Saves lots of distance calc for every Vs every in next step.



    // double overlaps = CalculateAtomicOverlaps(&enzyme, assembly);
    // std::cout << CalculateAtomicOverlaps(&enzyme, assembly) << std::endl;
    // std::ofstream myfile;

    //  std::stringstream ss1;
    //  ss1 << working_Directory + "/outputs/Overlaps.txt";
    //   std::cout << ss1.str() << std::endl;
    //   myfile.open (ss1.str(), std::ios::out | std::ios::app);
    //   std::string filename = assembly->GetName() + ": ";
    //   myfile << filename;
    //   myfile << overlaps << "\n";
    //   myfile.close();

    //CalculateAtomicOverlapsForAssemblyVector(&enzyme, substrate_snapshots);
    //double overlaps = enzyme->CalculateAtomicOverlaps(alsoMoving_atoms);
  //  std::cout << "enzyme overlap with substrate is " << overlaps << std::endl;
  //  return overlaps;
    //return enzyme->CalculateAtomicOverlaps(substrate);

//}


/*
double CalculateAtomicOverlapsForAssemblyVector(Assembly *assemblyA, AssemblyVector assemblies)
{
    std::vector<std::thread> th;

    //int nr_threads = assemblies->size();

    //Launch a group of threads
    for(AssemblyVector::iterator it = assemblies.begin(); it != assemblies.end(); ++it)
    {
        Assembly *current_assembly = *it;
        th.push_back(std::thread(CalculateAtomicOverlaps, assemblyA, current_assembly));
    }

//    for (int i = 0; i < nr_threads; ++i)
    //{
    //    th.push_back(std::thread(CalculateAtomicOverlaps, assemblyA, assemblies ));
    //}

    //Join the threads with the main thread
    for(auto &t : th)
    {
        t.join();
    }
}
*/


