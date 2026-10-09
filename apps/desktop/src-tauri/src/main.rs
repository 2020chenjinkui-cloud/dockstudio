//! Tauri 壳：窗口容器 + 后端 sidecar 进程编排。
//! 职责：起 PyInstaller 冻结的 dockstudio-backend（监听 127.0.0.1:8000），
//! 前端仍走 http://127.0.0.1:8000/api；应用退出时杀掉子进程。
//! 前置：需要 Rust 工具链 + `npm i -D @tauri-apps/cli@^2`（M0 本机未装，CI 上装）。

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use tauri::Manager;
use tauri_plugin_shell::process::{CommandChild, CommandEvent};
use tauri_plugin_shell::ShellExt;

struct BackendChild(std::sync::Mutex<Option<CommandChild>>);

fn main() {
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .setup(|app| {
            let shell = app.shell();
            let sidecar = shell
                .sidecar("dockstudio-backend")
                .expect("backend sidecar not bundled");
            let (mut rx, child) = sidecar.spawn().expect("failed to spawn backend");
            app.manage(BackendChild(std::sync::Mutex::new(Some(child))));
            tauri::async_runtime::spawn(async move {
                while let Some(event) = rx.recv().await {
                    if let CommandEvent::Stdout(line) | CommandEvent::Stderr(line) = event {
                        println!("[backend] {}", String::from_utf8_lossy(&line));
                    }
                }
            });
            Ok(())
        })
        .on_window_event(|window, event| {
            if let tauri::WindowEvent::Destroyed = event {
                let child = {
                    let state = window.state::<BackendChild>();
                    let x = state.0.lock().unwrap().take();
                    x
                };
                if let Some(child) = child {
                    let _ = child.kill();
                }
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
