pythonExec='/Applications/ParaView-5.5.1.app/Contents/bin/pvpython'
pythonScript='readlinesource_2D.py'
$pythonExec $pythonScript horizontal_line_${refLevName}.e
$pythonExec $pythonScript vertical_line_${refLevName}.e
