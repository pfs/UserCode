from hepdata_lib import Submission
submission = Submission()

from hepdata_lib import Table
table = Table("pa all")
table.description = "description."
table.location = "upper left."
table.keywords["observables"] = ["pa"]

from hepdata_lib import RootFileReader
reader = RootFileReader("root://eosuser.cern.ch//eos/user/v/vveckaln/analysis_MC13TeV_TTJets/plots/plotter.root")
Data = reader.read_hist_1d("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0/L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0")
Unc = reader.read_hist_1d("L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0/L_pull_angle_allconst_reco_leading_jet_scnd_leading_jet_DeltaRgt1p0_totalMCUncShape")

from hepdata_lib import Variable, Uncertainty

mmed = Variable("pa", is_independent=True, is_binned=False, units="rad")
mmed.values = signal["x"]

data = Variable("N", is_independent=False, is_binned=False, units="")
data.values = Data["y"]

unc = Uncertainty("Total", is_symmetric=True)
unc.values = Unc["dy"]
data.add_uncertainty(unc)

table.add_variable(mmed)
table.add_variable(data)

submission.add_table(table)

submission.create_files("example_output")

