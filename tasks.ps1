# よく使う処理をまとめたスクリプト（Makefile の代わり）
#
#   .\tasks.ps1 server                       結果サイトを作ってローカルサーバーで開く（http://127.0.0.1:8000/）
#   .\tasks.ps1 site                         結果サイトを作るだけ
#   .\tasks.ps1 benchmark                    合成データで全アプローチを比較
#   .\tasks.ps1 analyze data\hemibrain_mb.pkl
#                                            全アプローチ → 一致ペアの集計 → タイプの組の集計 → サイト作成
#
# 環境変数で設定を変えられる:  $env:PORT = 8080;  $env:N3D = 50
# 「スクリプトの実行が無効」と言われたら:
#   powershell -ExecutionPolicy Bypass -File .\tasks.ps1 server
param(
    [Parameter(Position = 0)][string]$Task = "help",
    [Parameter(Position = 1)][string]$Data = "data\hemibrain_mb.pkl"
)
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot
$py = Join-Path $PSScriptRoot ".venv\Scripts\python.exe"
$port = if ($env:PORT) { $env:PORT } else { "8000" }
$n3d = if ($env:N3D) { $env:N3D } else { "30" }

function Run([string[]]$cmd) {
    Write-Host "> python $($cmd -join ' ')" -ForegroundColor Cyan
    & $py @cmd
    if ($LASTEXITCODE -ne 0) { throw "失敗しました: $($cmd -join ' ')" }
}

switch ($Task) {
    "server"    { Run @("serve.py", "--port", $port, "--n3d", $n3d) }
    "site"      { Run @("build_site.py", "--n3d", $n3d) }
    "benchmark" { Run @("run_benchmark.py") }
    "analyze" {
        $stem = [IO.Path]::GetFileNameWithoutExtension($Data)
        Run @("run_all.py", "--data", $Data, "--compartments", "navis")
        Run @("compare_methods.py", "--results", "results\$stem", "--data", $Data, "--top-pct", "1", "--min-methods", "4")
        Run @("type_pairs.py", "--data", $Data)
        Run @("build_site.py", "--n3d", $n3d)
    }
    default {
        Get-Content $PSCommandPath | Select-Object -Skip 1 -First 10 | ForEach-Object { $_ -replace '^# ?', '' }
    }
}
