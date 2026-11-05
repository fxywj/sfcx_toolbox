' ============================================================
'  SFCX web GUI - double-click launcher (no console window)
'  Double-click me: starts the server silently and opens the
'  browser.  If it is already running, just opens the browser.
' ============================================================
Option Explicit
Dim sh, fso, base, port, logPath
port = 8765
Set sh = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
base = fso.GetParentFolderName(WScript.ScriptFullName)
logPath = fso.BuildPath(base, "sfcx_web_log.txt")

If ServerUp() Then
  sh.Run "http://127.0.0.1:" & port & "/", 1, False
  WScript.Quit 0
End If

Dim py
py = PickPython()
If py = "" Then
  MsgBox "No Python with numpy / scipy / matplotlib was found." & vbCrLf & vbCrLf & _
         "Install Python 3.9+ with those packages, or edit the" & vbCrLf & _
         "candidate list at the top of this file.", 16, "SFCX"
  WScript.Quit 1
End If

Dim cmd
cmd = "%COMSPEC% /c cd /d """ & base & """" & _
      " && set ""PYTHONPATH=" & base & """" & _
      " && """ & py & """ -m sfcx.webapp --port " & port & _
      " --open > """ & logPath & """ 2>&1"
sh.Run cmd, 0, False

Dim i
For i = 1 To 50
  WScript.Sleep 500
  If ServerUp() Then WScript.Quit 0
Next

Dim t
t = ""
On Error Resume Next
t = fso.OpenTextFile(logPath, 1).ReadAll
On Error GoTo 0
If Len(t) > 1200 Then t = Right(t, 1200)
MsgBox "The SFCX server did not start within 25 s." & vbCrLf & vbCrLf & t, 16, "SFCX"

Function ServerUp()
  ServerUp = False
  On Error Resume Next
  Dim h
  Set h = CreateObject("MSXML2.XMLHTTP")
  h.open "GET", "http://127.0.0.1:" & port & "/api/catalog", False
  h.send
  If h.status = 200 Then ServerUp = True
  On Error GoTo 0
End Function

Function PickPython()
  PickPython = ""
  Dim cands, c, e, waited
  ' add your own interpreter path to this list if needed
  cands = Array("python", "py")
  For Each c In cands
    On Error Resume Next
    Set e = sh.Exec("""" & c & """ -c ""import numpy,scipy,matplotlib""")
    If Err.Number = 0 Then
      waited = 0
      Do While e.Status = 0 And waited < 15000
        WScript.Sleep 100
        waited = waited + 100
      Loop
      If e.Status <> 0 Then
        If e.ExitCode = 0 Then
          On Error GoTo 0
          PickPython = c
          Exit Function
        End If
      End If
    End If
    Err.Clear
    On Error GoTo 0
  Next
End Function
