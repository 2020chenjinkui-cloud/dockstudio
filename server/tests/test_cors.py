"""CORS 回归：Tauri 生产 webview 来源必须放行，否则桌面端健康检查被浏览器拦截一直转圈。"""

from __future__ import annotations

ALLOWED = [
    "http://localhost:5173",
    "http://127.0.0.1:5173",
    "http://tauri.localhost",
    "tauri://localhost",
]


def test_tauri_webview_origins_allowed(client) -> None:
    for origin in ALLOWED:
        r = client.options(
            "/api/health",
            headers={"Origin": origin, "Access-Control-Request-Method": "GET"},
        )
        assert r.status_code == 200, f"{origin}: {r.status_code} {r.text}"
        assert r.headers.get("access-control-allow-origin") == origin


def test_unknown_origin_rejected(client) -> None:
    r = client.options(
        "/api/health",
        headers={"Origin": "https://evil.example", "Access-Control-Request-Method": "GET"},
    )
    assert r.status_code == 400
