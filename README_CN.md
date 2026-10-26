# 再入黑障状态估计

[English](README.md) · [任务指标](tasks.md)

Rust导航估计研究实现：24维误差状态平方根EKF、因果角速度观测器、按测量时刻
回退重放，以及固定时滞RTS平滑。支持GNSS位置/速度解和原始伪距/多普勒观测。

## 参考场景

命令行默认 `SCENARIO=reference-entry`，使用显式合成的下降轨迹与气动表；原有
气动、RCS预览数据不足以标定完整飞行器模型。`SCENARIO=legacy` 保留原轨迹，
库调用的 `SimConfig::default()` 保留旧初始条件。

IMU/GNSS/星敏/雷达频率分别为400/10/2/20 Hz。星敏和雷达分别延迟120/70 ms，
雷达仅在25 km以下有效。GNSS在 `[120,160) s` 中断；星敏按测量时刻在
`[130,150] s` 不可用。默认含5%重尾污染和20%随机缺测。
真值噪声与滤波器调参分离；同种子切换模式保持物理真值和IMU随机流一致。

## 构建与运行

```bash
cargo test --release --locked
cargo test --release --locked --no-default-features --features flight
cargo build --release --locked --bin harness --bin repro --bin resourcebench
ASSESSMENT=causal TIGHT=1 OUT_DIR=output_tight target/release/harness 1000 --check
ASSESSMENT=deadline OUT_DIR=output_pvt target/release/harness 1000 --check
TIGHT=1 target/release/repro 0
target/release/resourcebench --tight
```

`--check` 在所选数值验收失败时返回非零状态。默认 `ASSESSMENT=causal` 检查
160–162秒的801个即时输出；`ASSESSMENT=deadline` 只使用截至162秒已收到的
观测，平滑同一段的201个检查点。两种结果始终分开保存，平滑结果不能冒充即时输出。
`SEED_OFFSET` 指定批量试验的起始种子。

`scripts/validate-reference.sh` 运行测试、两模式验收和资源探针。CSV、图表和
运行日志保存在已被Git忽略的 `output*` 目录中。

[设计说明](docs/design.md)、[接口约定](docs/interface.md)和
[验证流程](docs/verification.md)定义实现与使用方法。
资源探针分别报告常规周期和RTS批处理耗时；在ARM参考目标上使用 `--check-time`
检查5 ms墙钟预算。合成模型与主机测量不能代替真实飞行器和目标硬件验证。
