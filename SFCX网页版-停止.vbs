' ============================================================
'  SFCX web GUI - stop the hidden server
'  Double-click me to close the background SFCX server.
' ============================================================
Option Explicit
Dim sh, port, portStr
port = 8765
Set sh = CreateObject("WScript.Shell")
portStr = ":" & port

Dim ex, line, toks, pid, killed
Set ex = sh.Exec("%COMSPEC% /c netstat -ano -p tcp | findstr " & portStr)
Do While ex.Status = 0
  WScript.Sleep 50
Loop
killed = 0
Do While Not ex.StdOut.AtEndOfStream
  line = Trim(ex.StdOut.ReadLine)
  If InStr(line, "127.0.0.1:" & port) > 0 And InStr(line, "LISTENING") > 0 Then
    toks = Split(line)
    pid = toks(UBound(toks))
    If IsNumeric(pid) Then
      If CInt(pid) <> 0 Then
        sh.Run "%COMSPEC% /c taskkill /F /PID " & pid, 0, True
        killed = killed + 1
      End If
    End If
  End If
Loop

If killed > 0 Then
  MsgBox "SFCX server stopped (" & killed & " process).", 64, "SFCX"
Else
  MsgBox "No SFCX server is running on port " & port & ".", 48, "SFCX"
End If
