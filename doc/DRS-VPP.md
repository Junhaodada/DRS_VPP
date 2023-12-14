## 标题

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701433308472-a5353445-2a22-436b-bd03-476a510370bd.png#averageHue=%23e5e2e0&clientId=u328af614-f285-4&from=paste&height=160&id=udb4fdb27&originHeight=160&originWidth=887&originalType=binary&ratio=1&rotation=0&showTitle=false&size=54920&status=done&style=none&taskId=u30f541b8-f51d-4826-bacf-bdc4f47583d&title=&width=887)

## 摘要

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701433744427-77bee0f4-2dcd-4bf4-885c-605644c16393.png#averageHue=%23f4f1ee&clientId=u328af614-f285-4&from=paste&height=543&id=ua425811c&originHeight=543&originWidth=1005&originalType=binary&ratio=1&rotation=0&showTitle=false&size=221659&status=done&style=none&taskId=u412624d2-1269-4677-9a82-5253720b254&title=&width=1005)

研究背景

* 化石燃料对环境的污染，需要面向不同的用户级别，设计综合能源系统，向用户供给电、冷能等，系统的设计目标是保护环境+节约成本。

研究目标

* 分布式可再生能源在PIES的高效管理

贡献

* 考虑PV和负载不确定对VPP调度的影响，提出两阶段鲁棒优化模型

## 引言

## VPP框架

### VPP物理模型

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701435299654-377abba2-5283-4cec-8945-abf94e58e807.png#averageHue=%23f6f4f1&clientId=u328af614-f285-4&from=paste&height=692&id=u47802ef5&originHeight=692&originWidth=1032&originalType=binary&ratio=1&rotation=0&showTitle=false&size=559070&status=done&style=none&taskId=u07fd45c3-71db-4170-bfe8-55f4f89ca63&title=&width=1032)

能量

* 电能
   - 来源：PV, ESS, the micro-turbine generator (MTG) and the main grid
* 冷能
   - 来源： electric chillers (EC)  

研究重点

* ** two sources of uncertainty：the PVs output** and **the load demands**. 
* **Reducing the impact of uncertainties on scheduling results** is the focus of the differential modeling of uncertainties

PV输出不确定性

* 原因：由于天气情况导致的不确定性，光照辐射和云量
* 措施：probability density function (PDF) of PV to estimate

负载不确定性

* 对于电负荷，我们考虑根据综合需求响应(IDR)调整用电量，以进一步减轻高发电量对电价、运营成本的不利影响。
* 冷负荷，它通常基于建筑物的温度和热力学模型进行建模

改进

* DL改进不确定性
* 碳排放加惩罚

### VPP数学模型

需求响应

1. 电负载模型
* 两种负载类型：transferable loads (TLs) and interruptible loads (ILs)
* $P^{TL}_t$代表用电转移

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701516259817-54da611c-5c21-4214-a6a5-c78a0cc0b406.png#averageHue=%23fbfaf9&clientId=u4a00da34-0014-4&from=paste&height=217&id=ud3c1bf43&originHeight=325&originWidth=820&originalType=binary&ratio=1.5&rotation=0&showTitle=false&size=51956&status=done&style=none&taskId=u72348dce-8445-4f6f-9373-cb135179d5b&title=&width=546.6666666666666)

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701516605219-f4fb0eed-4e32-41ed-b2b8-4e1b6a05c844.png#averageHue=%23fcfbfa&clientId=u4a00da34-0014-4&from=paste&height=57&id=u924c5d6d&originHeight=85&originWidth=877&originalType=binary&ratio=1.5&rotation=0&showTitle=false&size=11703&status=done&style=none&taskId=ua9f74366-440d-4a4b-a1ec-20d7e4f52cd&title=&width=584.6666666666666)

2. 冷负载

冷负载分为：工厂冷负载、居民区和办公楼冷负载两种模型。

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701517548973-55279b79-0352-45c1-b096-de58fe4579ff.png#averageHue=%23f2efec&clientId=u4a00da34-0014-4&from=paste&height=29&id=uf56fe69b&originHeight=43&originWidth=730&originalType=binary&ratio=1.5&rotation=0&showTitle=false&size=15831&status=done&style=none&taskId=u6d7e2577-1b63-450f-8f3f-1f1e703936c&title=&width=486.6666666666667)

工厂冷负载

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701517485612-4c98f03e-2778-4429-9075-f7acb2518d9e.png#averageHue=%23fbf9f8&clientId=u4a00da34-0014-4&from=paste&height=110&id=ue7fc9619&originHeight=165&originWidth=969&originalType=binary&ratio=1.5&rotation=0&showTitle=false&size=37969&status=done&style=none&taskId=ucdebd446-06c4-4fd7-b863-463d5b9a1cd&title=&width=646)

居民区和办公楼冷负载

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701517593932-e89c7663-9bee-47ae-a207-9259f5f49f1d.png#averageHue=%23faf9f8&clientId=u4a00da34-0014-4&from=paste&height=205&id=ud2a4c8f6&originHeight=307&originWidth=990&originalType=binary&ratio=1.5&rotation=0&showTitle=false&size=71229&status=done&style=none&taskId=u4eeb21f0-d273-4e1e-aa6a-cb7e92ef810&title=&width=660)

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701517830668-f7d19e84-96b3-4a4b-9a94-7287018f3d4c.png#averageHue=%23fbfaf9&clientId=u4a00da34-0014-4&from=paste&height=53&id=u8758828f&originHeight=79&originWidth=865&originalType=binary&ratio=1.5&rotation=0&showTitle=false&size=17120&status=done&style=none&taskId=u7191bdee-ad9d-47a6-8a3e-69e8e191e94&title=&width=576.6666666666666)

![image.png](https://cdn.nlark.com/yuque/0/2023/png/1709231/1701517860599-f705fcea-5a88-4080-bcb4-8f8927866bf0.png#averageHue=%23fbfafa&clientId=u4a00da34-0014-4&from=paste&height=181&id=uf10f70c8&originHeight=271&originWidth=873&originalType=binary&ratio=1.5&rotation=0&showTitle=false&size=47505&status=done&style=none&taskId=u7f49efd7-fbd3-4e5b-8bc1-d95a8f8ebfe&title=&width=582)

成本
生产维护成本
发电成本

## GAN

输出N组PV输出和N组冷负载场景

## VPP优化模型

![](https://cdn.nlark.com/yuque/0/2023/jpeg/1709231/1702179509893-d7991ac9-d39b-491e-b813-a2ec6a44dabd.jpeg)

## 模型转换和求解

## 数值实验

## 知识点

[燃气涡轮发电机](https://baike.baidu.com/item/%E7%87%83%E6%B0%94%E6%B6%A1%E8%BD%AE%E5%8F%91%E7%94%B5%E6%9C%BA/5938147)
