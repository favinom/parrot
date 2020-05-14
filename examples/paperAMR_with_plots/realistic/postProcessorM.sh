pythonExec='/Applications/ParaView-5.5.1.app/Contents/bin/pvpython'

pythonScript='readlinesource_2D.py'

$runline plotLine_master_unresolved.i mRefType=${refinementType} mRefLev=${refLev}
# $pythonExec $pythonScript horizontal_line_${refinementType}_${refLev}.e
# $pythonExec $pythonScript vertical_line_${refinementType}_${refLev}.e

