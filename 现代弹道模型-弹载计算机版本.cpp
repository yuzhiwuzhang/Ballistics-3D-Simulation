/* ============================================================
   现代弹道模型 - 完整独立C语言版本
   编译：gcc ballistics.c -o ballistics.exe -lm
   运行：ballistics.exe
   ============================================================ */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* -------------------- 物理常数 -------------------- */
#define RHO0        1.225f
#define T0          288.15f
#define P0          101325.0f
#define L           0.0065f
#define R           287.05f
#define G0          9.80665f
#define RE          6371000.0f
#define OMEGA       0.00007292f
#define GAMMA       1.4f

/* -------------------- 大气查找表（0-10000米，步长100米） -------------------- */
#define ATMO_SIZE   101

typedef struct {
    float rho;
    float a;
    float T;
} AtmoData;

/* 只保留关键数据点（实际使用时需补全101行，这里用公式替代） */
static float get_atmo_rho(float h) {
    /* 用标准大气公式替代查表，这样就不用写101行数据了 */
    if (h < 0) h = 0;
    float T = T0 - L * h;
    float P = P0 * powf(T / T0, G0 / (R * L));
    return P / (R * T);
}

static float get_atmo_a(float h) {
    if (h < 0) h = 0;
    float T = T0 - L * h;
    return sqrtf(GAMMA * R * T);
}

/* -------------------- 阻力系数表（0-3马赫，步长0.1） -------------------- */
static const float cd_table[] = {
    0.120f, 0.120f, 0.121f, 0.123f, 0.126f, 0.131f, 0.138f,
    0.148f, 0.162f, 0.180f, 0.220f, 0.290f, 0.420f, 0.550f,
    0.600f, 0.620f, 0.630f, 0.640f, 0.645f, 0.650f, 0.655f,
    0.658f, 0.660f, 0.662f, 0.664f, 0.666f, 0.668f, 0.670f,
    0.672f, 0.674f, 0.676f, 0.678f
};
#define CD_TABLE_SIZE (sizeof(cd_table) / sizeof(cd_table[0]))

static float get_cd(float mach) {
    if (mach <= 0.0f) return cd_table[0];
    if (mach >= 3.0f) return cd_table[CD_TABLE_SIZE - 1];
    int idx = (int)(mach / 0.1f);
    float frac = (mach - idx * 0.1f) / 0.1f;
    if (idx >= CD_TABLE_SIZE - 1) return cd_table[CD_TABLE_SIZE - 1];
    return cd_table[idx] + frac * (cd_table[idx + 1] - cd_table[idx]);
}

/* -------------------- 弹道计算核心 -------------------- */
void compute_trajectory(
    float v0, float theta_deg, float lat_deg, float drag_coeff,
    float max_time, float dt,
    float *range, float *max_h, float *flight_time,
    int *steps
) {
    float theta = theta_deg * 3.14159265f / 180.0f;
    float lat = lat_deg * 3.14159265f / 180.0f;
    float sin_lat = sinf(lat);

    float x = 0.0f, y = 0.0f;
    float vx = v0 * cosf(theta);
    float vy = v0 * sinf(theta);

    *max_h = 0.0f;
    *range = 0.0f;
    *flight_time = 0.0f;
    *steps = 0;

    float t = 0.0f;

    while (t < max_time) {
        /* 1. 大气参数 */
        float rho = get_atmo_rho(y);
        float a = get_atmo_a(y);

        /* 2. 速度 */
        float v = sqrtf(vx * vx + vy * vy);
        if (v < 0.001f) break;

        /* 3. 马赫数 + 阻力系数 */
        float mach = v / a;
        float Cd = get_cd(mach);

        /* 4. 有效阻力 */
        float drag = drag_coeff * (rho / RHO0) * (Cd / 0.15f);

        /* 5. 加速度 */
        float g_h = G0 * (RE / (RE + y)) * (RE / (RE + y));
        float ax = -drag * vx * v + 2.0f * OMEGA * vy * sin_lat;
        float ay = -drag * vy * v - g_h - 2.0f * OMEGA * vx * sin_lat;

        /* 6. 改进欧拉积分（中点法） */
        float vx_mid = vx + ax * dt * 0.5f;
        float vy_mid = vy + ay * dt * 0.5f;

        float v_mid = sqrtf(vx_mid * vx_mid + vy_mid * vy_mid);
        float mach_mid = v_mid / a;
        float Cd_mid = get_cd(mach_mid);
        float drag_mid = drag_coeff * (rho / RHO0) * (Cd_mid / 0.15f);

        float ax_mid = -drag_mid * vx_mid * v_mid + 2.0f * OMEGA * vy_mid * sin_lat;
        float ay_mid = -drag_mid * vy_mid * v_mid - g_h - 2.0f * OMEGA * vx_mid * sin_lat;

        vx += ax_mid * dt;
        vy += ay_mid * dt;
        x += vx_mid * dt;
        y += vy_mid * dt;
        t += dt;

        /* 记录 */
        (*steps)++;
        if (y > *max_h) *max_h = y;

        /* 落地判断 */
        if (y < 0.0f) {
            *range = x;
            *flight_time = t;
            return;
        }
    }

    /* 超出最大时间 */
    *range = x;
    *flight_time = t;
}

/* -------------------- 主程序（控制台输出） -------------------- */
int main() {
    printf("============================================================\n");
    printf("         现代弹道模型 - C语言弹载级实现\n");
    printf("         编译时间: 2026年8月\n");
    printf("============================================================\n\n");

    /* 参数设置（模拟155mm榴弹炮） */
    float v0 = 700.0f;          /* 初速度 m/s */
    float theta = 45.0f;        /* 发射角 度 */
    float lat = 34.0f;          /* 纬度（邳州/徐州） */
    float drag_coeff = 0.0005f; /* 阻力系数 */
    float max_time = 40.0f;     /* 最大时间 s */
    float dt = 0.01f;           /* 时间步长 s */

    float range, max_h, flight_time;
    int steps;

    printf("输入参数:\n");
    printf("  初速度: %.1f m/s\n", v0);
    printf("  发射角: %.1f 度\n", theta);
    printf("  纬度: %.1f 度\n", lat);
    printf("  阻力系数: %.6f\n", drag_coeff);
    printf("  时间步长: %.3f s\n\n", dt);

    compute_trajectory(v0, theta, lat, drag_coeff, max_time, dt,
                       &range, &max_h, &flight_time, &steps);

    printf("计算结果:\n");
    printf("  射程: %.1f 米 (%.2f 公里)\n", range, range / 1000.0f);
    printf("  最大高度: %.1f 米\n", max_h);
    printf("  飞行时间: %.2f 秒\n", flight_time);
    printf("  计算步数: %d\n", steps);

    /* ---- 对比：不同角度下的射程 ---- */
    printf("\n============================================================\n");
    printf("角度-射程关系表 (初速度 %.0f m/s):\n", v0);
    printf("============================================================\n");
    printf("  角度(°) |  射程(km) |  最大高度(m) |  飞行时间(s)\n");
    printf("------------------------------------------------------------\n");

    float angles[] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65};
    int num_angles = sizeof(angles) / sizeof(angles[0]);

    for (int i = 0; i < num_angles; i++) {
        float th = angles[i];
        float rng, hmax, tfly;
        int stp;
        compute_trajectory(v0, th, lat, drag_coeff, 60.0f, dt,
                           &rng, &hmax, &tfly, &stp);
        printf("  %6.0f   |  %7.2f   |  %8.0f     |  %7.2f\n",
               th, rng / 1000.0f, hmax, tfly);
    }

    printf("============================================================\n");
    printf("说明: 此代码使用改进欧拉积分 + 标准大气模型\n");
    printf("      阻力系数随马赫数变化（跨音速区域已建模）\n");
    printf("      含科里奥利力修正（地球自转）\n");
    printf("============================================================\n");

    return 0;
}
