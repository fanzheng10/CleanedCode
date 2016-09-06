use NAMD;
use GENERAL;
use DEFINITIONS;
use Getopt::Long;
use File::Spec;

my %opts = getInput();

# file naming conventions
my $ipdbf = "_init.pdb";
my $opdbf = "_start.pdb";
my $opsff = GENERAL::GetBase($opdbf) . ".psf";
my $wpdbf = GENERAL::GetBase($opdbf) . "_wat.pdb";
my $wpsff = GENERAL::GetBase($wpdbf) . ".psf";
my $pdbf = $opdbf; my $psff = $opsff;

# create local space
my $odir = GENERAL::GetDir();
my $ldir = GENERAL::createLocalSpace($opts{tmp});
GENERAL::csystem("cp $opts{p} $ldir/$ipdbf");
GENERAL::cchdir($ldir);

# create original PSF
NAMD::createPSF($ipdbf, $opts{top}, $opdbf, $opsff, 'genopts' => $opts{psfopts}, 'patches' => $opts{patches}, 'vmd_bin' => '/home/grigoryanlab/local/bin/vmd');

# solvate
my $S;
if (!defined($opts{is})) {
  $S = NAMD::waterBoxVMD($opdbf, $opsff, $opts{pad}, $wpdbf, $wpsff);
  $pdbf = $wpdbf; $psff = $wpsff;
  $S->{xsize} *= 0.95; $S->{ysize} *= 0.95; $S->{zsize} *= 0.95; # because the box will compress (this adjustment should help avoid the "FATAL ERROR: Periodic cell has become too small for original patch grid!"
  NAMD::addIons($pdbf, $psff, GENERAL::GetBase($pdbf), "neutralize");
}

# start preparing NAMD script
my $bname = "namdrun_run";
my $namd = NAMD::new('exec' => "/home/grigoryanlab/library/bin/namd2", 'cmdfile' => "$bname.ctl", 'keepout' => 1, 'np' => $opts{np}, 'charmrun' => "/home/grigoryanlab/library/bin/charmrun");
$namd->setParameter("forcefield", "parameters", $opts{par});
$namd->setParameter("forcefield", "cutoff", $opts{cut});
$namd->setParameter("forcefield", "switching", "on");
$namd->setParameter("forcefield", "switchdist", $opts{switch});
$namd->setParameter("forcefield", "pairlistdist", $opts{cut} + 1.5);
$namd->setParameter("variables", "set temperature", $opts{t}, "temperature in Kelvin", "default_f", 1);
$namd->prependParameterGroup("structure");
$namd->appendParameter("structure", "structure", $psff);
$namd->appendParameter("structure", "coordinates", $pdbf);
$namd->appendParameter("integrator", "seed", $opts{seed}, "random number generator seed") if (defined($opts{seed}));
$namd->setParameter("integrator", "timestep", $opts{ts});

# solvent stuff
if (!defined($opts{is})) {
  $namd->appendParameterGroup("boundary conditions");
  $namd->appendParameter("boundary conditions", "cellBasisVector1", "$S->{xsize} 0 0");
  $namd->appendParameter("boundary conditions", "cellBasisVector2", "0 $S->{ysize} 0");
  $namd->appendParameter("boundary conditions", "cellBasisVector3", "0 0 $S->{zsize}");
  $namd->appendParameter("boundary conditions", "cellOrigin", "$S->{cx} $S->{cy} $S->{cz}");
  $namd->appendParameter("boundary conditions", "wrapAll", "on");

  if (defined($opts{Pext})) {
    $namd->appendParameterGroup("constant pressure");
    $namd->appendParameter("constant pressure", "useGroupPressure", "yes");
    $namd->appendParameter("constant pressure", "useFlexibleCell", "no");
    $namd->appendParameter("constant pressure", "useConstantRatio", "no");
    $namd->appendParameter("constant pressure", "useConstantArea", "no");
    $namd->appendParameter("constant pressure", "langevinPiston", "on");
    $namd->appendParameter("constant pressure", "langevinPistonTarget", $opts{Pext}, "pressure in bar");
    $namd->appendParameter("constant pressure", "langevinPistonPeriod", "100");
    $namd->appendParameter("constant pressure", "langevinPistonDecay", "50");
    $namd->appendParameter("constant pressure", "langevinPistonTemp", "\$temperature");
  }

  if (defined($opts{pme})) {
    $namd->appendParameter("integrator", "PME", "on");
    $namd->appendParameter("integrator", "PMEInterpOrder", 4, "NAMD manual default");
    $namd->appendParameter("integrator", "PMEGridSpacing", 1);
  }
} else {
  $namd->appendParameterGroup("implicit solvent");
  $namd->appendParameter("implicit solvent", "GBIS", "on");
  $namd->appendParameter("implicit solvent", "SASA", "on");
}

# fixed atoms, if any
if (defined($opts{f})) {
  my $pdb = PDB::new($pdbf);
  my @atoms = $pdb->conAtoms();
  my $n = 0;
  foreach my $atom (@atoms) {
    if (PDB::atomStr($atom) =~ /$opts{f}/) {
      $atom->{B} = 1;
      $n++;
    } else {
      $atom->{B} = 0;
    }
  }
  my $fpdbf = "_fixed.pdb";
  $pdb->writePDB($fpdbf, "");
  $namd->appendParameterGroup("fixed atoms");
  $namd->appendParameter("fixed atoms", "fixedAtoms", "on");
  $namd->appendParameter("fixed atoms", "fixedAtomsForces", "off");
  $namd->appendParameter("fixed atoms", "fixedAtomsFile", $fpdbf);
  $namd->appendParameter("fixed atoms", "fixedAtomsCol", "B");
  printf("Fixing %d atoms...\n", $n);
}

# IO
$namd->setParameter("io", "outputEnergies", $opts{ns});
$namd->appendParameter("io", "dcdfreq", $opts{ns});

# commands
$namd->command("minimize $opts{mini}") if (defined($opts{mini}));
$namd->command("run $opts{n}");

# run NAMD
if (!defined($opts{dry})) {
  my $ret = $namd->runNAMD();
  printf("Return value from first run '$ret'\n");
} else {
  $namd->writeCTLFile();
}
GENERAL::cchdir($odir);

if (defined($opts{odir})) {
  GENERAL::cmkdir($opts{odir}) if (!-e $opts{odir});
  GENERAL::csystem("cp -r $ldir/* $opts{odir}/");
} else {
  GENERAL::csystem("cp $ldir/$bname.namdout.dcd $opts{ob}.dcd") if (-e "$ldir/$bname.namdout.dcd");
  GENERAL::csystem("cp $ldir/$pdbf $opts{ob}.pdb") if (-e "$ldir/$pdbf");
  GENERAL::csystem("cp $ldir/$psff $opts{ob}.psf") if (-e "$ldir/$psff");
}

# clean up
GENERAL::destroyLocalSpace($ldir);

sub getInput {
  my $usage = GENERAL::usage(
            "Performs an MD simulation on a PDB structure with NAMD.", "Required switches:", "",
            "-p", "PDB file name",
            "-n", "number of steps of MD (in fs)",
            "Optional switches:", "",
            "--f", "selecton of atoms to keep fixed. Regular expression applied to atom strings (format: <chain>_<residue name>_<residue number>_<atom name>).",
            "--ns", "sampling frequency: number of steps between .dcd writing/energy printing. Default is 100.",
            "--ts", "integration timestep in fs. Default is 1 fs.",
            "-t", "temperature for the simulation. The default is 298.15 K.",
            "--is", "implicit solvent (do not solvate with waters, use GBSA)",
            "--pad", "for explicit-solvent simulations -- the padding (in Angstroms) for how much water to add from each side (default is 5 A)",
            "--Pext", "constant pressure to keep (in bar; 1.01325 is 1 atm). If not specified, constant pressure will not be applied.",
            "--mini", "Pre-minimize for this many steps of steepest descent.",
            "--cut", "interaction cutoff distance (default is 12)",
            "--switch", "distance where switching function picks up (default is 10)",
            "--pme", "use Particle Mesh Ewald (PME) summation for long-range electrostatics (default is not to use PME). With PME, could use shorter cutoffs (e.g., --cut 10 and --switch 6).",
            "--seed", "initial integer seed for assigning velocities (by default, a random seed is used by combining a random integer and the PID)",
            "--top", "custom topology file",
            "--par", "custom parameter file",
            "--ob", "base name for various kinds of output",
            "--patches", "a list of patches to apply (i.e. 'patch DISU BPTI:5 BPTI:55 | patch NTER A:1' - one set of patches per chain, separated by |)",
            "--psfopts", "aditional options for generating PSF (e.g. 'first none; last none | first nter; last cter' - one set of options per chain, separated by |)",
            "--np", "number of processes to use (default is 1)",
            "--tmp", "temporary directory path prefix",
            "--odir", "path to copy the entire working directory to",
            "--dry", "make all of the preparations, but do not actually run NAMD");

  my %opts;
  my (@rType, @rOpts);
  GetOptions (\%opts, "p=s", "n=s", "ns=s", "t=s", "is", "mini=s", "seed=s", "top=s", "par=s", "ob=s", "patches=s", "tmp=s", "psfopts=s", "pad=s", "np=s", "Pext=s", "odir=s", "f=s", "dry", "ts=s", "pme", "cut=s", "switch=s");
  # Required options
  if (!(defined($opts{'p'}) && defined($opts{'n'}))) {
    die($usage);
  }
  GENERAL::assert((-e $opts{p}) ? 1 : 0, "file $opts{p} not found");
  $opts{top} = File::Spec->rel2abs($opts{top}) if (defined($opts{top}));
  $opts{par} = File::Spec->rel2abs($opts{par}) if (defined($opts{par}));
  $opts{top} = "/u/gevorg/work/ProCEDe/branches/gevorg-branch/param/toppar_c35b2_c36a2/top_all22_prot_namd.inp" if (!defined($opts{top}));
  $opts{par} = "/u/gevorg/work/ProCEDe/branches/gevorg-branch/param/toppar_c35b2_c36a2/par_all22_prot.inp" if (!defined($opts{par}));
  GENERAL::assert((-e $opts{top}) ? 1 : 0, "file $opts{top} not found");
  GENERAL::assert((-e $opts{par}) ? 1 : 0, "file $opts{par} not found");
  GENERAL::assert(GENERAL::isInteger($opts{n}), "-n should be integer!");
  $opts{ns} = 100 if (!defined($opts{ns}));
  GENERAL::assert(GENERAL::isInteger($opts{ns}), "-ns should be integer!");
  $opts{ts} = 1 if (!defined($opts{ts}));
  GENERAL::assert(GENERAL::isInteger($opts{ts}) && ($opts{ts} > 0), "-ts should be a positive integer!");
  $opts{seed} = $$ + int(rand()*1000000) if (!defined($opts{'s'}));
  GENERAL::assert(GENERAL::isInteger($opts{seed}), "--seed should be integer!");
  $opts{t} = 298.15 if (!defined($opts{t}));
  GENERAL::assert(GENERAL::isNumeric($opts{t}) && ($opts{t} > 0), "-t should be numeric and positive!");
  GENERAL::assert(GENERAL::isInteger($opts{mini}), "--mini must be integer!") if (defined($opts{mini}));
  if (!defined($opts{ob})) { $opts{ob} = GENERAL::GetBase($opts{p}) . "_$opts{n}_out"; }
  $opts{pad} = 5 if (!defined($opts{pad}));
  GENERAL::assert(GENERAL::isNumeric($opts{pad}) && ($opts{pad} > 0), "--pad must be numeric and positive!");
  $opts{np} = 1 if (!defined($opts{np}));
  GENERAL::assert(GENERAL::isNumeric($opts{np}) && ($opts{np} > 0), "--np must be numeric and positive!");
  GENERAL::assert(GENERAL::isNumeric($opts{Pext}) && ($opts{Pext} > 0), "--Pext must be numeric and positive!") if (defined($opts{Pext}));
  $opts{cut} = 12 if (!defined($opts{cut}));
  $opts{switch} = 10 if (!defined($opts{switch}));
  GENERAL::assert(GENERAL::isNumeric($opts{cut}) && ($opts{cut} > 0), "--cut must be numeric and positive!");
  GENERAL::assert(GENERAL::isNumeric($opts{switch}) && ($opts{switch} > 0) && ($opts{switch} <= $opts{cut}), "--switch must be numeric, positive, and less than the cutoff distance!");
  foreach my $key ("patches", "psfopts") {
    if (defined($opts{$key})) {
      my @key = split(/\|/, $opts{$key});
      $opts{$key} = \@key;
    } else {
      $opts{$key} = undef;
    }
  }
  $opts{tmp} = "/data/scratch" if (!defined($opts{tmp}));

  return %opts;
}
