"""
WGA-GP 生成光伏数据
"""
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import matplotlib.pyplot as plt


class Generator(nn.Module):
    """定义生成器"""

    def __init__(self):
        super(Generator, self).__init__()
        self.model = nn.Sequential(
            # 1D 卷积层 1
            nn.Conv1d(1, 32, kernel_size=3, stride=1),
            nn.ReLU(),
            # 批归一化层 1
            nn.BatchNorm1d(32, momentum=0.8),
            # 1D 卷积层 2
            nn.Conv1d(32, 64, kernel_size=3, stride=1),
            nn.ReLU(),
            # 批归一化层 2
            nn.BatchNorm1d(64, momentum=0.8),
            # 1D 卷积层 3
            nn.Conv1d(64, 1, kernel_size=3, stride=1),
            nn.Tanh()
        )

    def forward(self, x):
        return self.model(x)


class Discriminator(nn.Module):
    """定义判别器"""

    def __init__(self):
        super(Discriminator, self).__init__()
        self.f1 = nn.Sequential(
            nn.Conv1d(1, 32, kernel_size=3, stride=1),
            nn.LeakyReLU(),
        )
        self.f2 = nn.Sequential(
            # 批归一化层 1
            nn.BatchNorm1d(32, momentum=0.8),
            # 1D 卷积层 2
            nn.Conv1d(32, 64, kernel_size=3, stride=1),
            nn.LeakyReLU(),
        )
        self.f3 = nn.Sequential(
            # 批归一化层 2
            nn.BatchNorm1d(64, momentum=0.8),
            # 1D 卷积层 3
            nn.Conv1d(64, 16, kernel_size=3, stride=1),
            nn.LeakyReLU(),
        )
        self.f4 = nn.Sequential(
            # 批归一化层 3
            nn.BatchNorm1d(16, momentum=0.8),

        )
        self.f5 = nn.Sequential(
            # 全连接层
            nn.Linear(18, 1)

        )
        # self.model = nn.Sequential(
        #     # 1D 卷积层 1
        #     nn.Conv1d(1, 32, kernel_size=3, stride=1),
        #     nn.LeakyReLU(),
        #     # 批归一化层 1
        #     nn.BatchNorm1d(32, momentum=0.8),
        #     # 1D 卷积层 2
        #     nn.Conv1d(32, 64, kernel_size=3, stride=1),
        #     nn.LeakyReLU(),
        #     # 批归一化层 2
        #     nn.BatchNorm1d(64, momentum=0.8),
        #     # 1D 卷积层 3
        #     nn.Conv1d(64, 16, kernel_size=3, stride=1),
        #     nn.LeakyReLU(),
        #     # 批归一化层 3
        #     nn.BatchNorm1d(16, momentum=0.8),
        #     # 全连接层
        #     nn.Linear(16, 1)
        # )

    def forward(self, x):
        # input [64,1,24]
        x1 = self.f1(x)
        # print(x1.shape)
        x2 = self.f2(x1)
        # print(x2.shape)
        x3 = self.f3(x2)
        # print(x3.shape)
        x4 = self.f4(x3)
        # print(x4.shape)
        x5 = self.f5(x4)
        # print(x5.shape)
        return x5


# 定义生成器和判别器的实例
generator = Generator()
discriminator = Discriminator()

# 使用 GPU 进行训练（如果可用）
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
generator.to(device)
discriminator.to(device)

# 定义损失函数和优化器
criterion = nn.BCEWithLogitsLoss()
optimizer_G = optim.Adam(generator.parameters(), lr=0.0002, betas=(0.5, 0.999))
optimizer_D = optim.Adam(discriminator.parameters(),
                         lr=0.0002, betas=(0.5, 0.999))

# 定义梯度惩罚函数


def compute_gradient_penalty(discriminator, real_samples, fake_samples):
    # print(real_samples.shape)
    # print(fake_samples.shape)
    alpha = torch.rand(real_samples.size(0), 1, 1, device=device)
    # print(alpha.shape)
    interpolates = (alpha * real_samples + (1 - alpha)
                    * fake_samples).requires_grad_(True)
    d_interpolates = discriminator(interpolates)
    fake = torch.ones(d_interpolates.size(), device=device)
    gradients = torch.autograd.grad(outputs=d_interpolates, inputs=interpolates,
                                    grad_outputs=fake, create_graph=True, retain_graph=True, only_inputs=True)[0]
    gradients = gradients.view(gradients.size(0), -1)
    gradient_penalty = ((gradients.norm(2, dim=1) - 1) ** 2).mean()
    return gradient_penalty


# 数据集
# 生成服从正态分布的随机数据作为真实数据样本
mean = 0  # 正态分布的均值
std = 1   # 正态分布的标准差
sample_size = 1000  # 样本数量
data_dim = 24  # 生成的数据维度
real_data = np.random.normal(mean, std, (sample_size, data_dim))

batch_size = 64
data_loader = DataLoader(TensorDataset(
    torch.Tensor(real_data)), batch_size=batch_size, shuffle=True, drop_last=True)

# 训练参数
n_epochs = 100
clip_value = 0.01
lambda_gp = 10  # 梯度惩罚系数

# 训练循环
G_losses, D_losses = [], []
for epoch in range(n_epochs):
    for batch in data_loader:
        real_samples = batch[0].to(device)
        fake_samples = generator(torch.randn(
            batch_size, 1, 30).to(device))

        # 计算判别器损失
        optimizer_D.zero_grad()
        real_samples = real_samples.unsqueeze(1)
        real_output = discriminator(real_samples)
        fake_output = discriminator(fake_samples.detach())
        loss_real = -real_output.mean()
        loss_fake = fake_output.mean()
        gradient_penalty = compute_gradient_penalty(
            discriminator, real_samples, fake_samples)
        loss_D = loss_real + loss_fake + lambda_gp * gradient_penalty
        loss_D.backward(retain_graph=True)
        optimizer_D.step()

        # 截断参数以防止梯度爆炸
        for p in discriminator.parameters():
            p.data.clamp_(-clip_value, clip_value)

        # 计算生成器损失
        optimizer_G.zero_grad()
        fake_output = discriminator(fake_samples)
        loss_G = -fake_output.mean()
        loss_G.backward(retain_graph=True)
        optimizer_G.step()

        G_losses.append(loss_G.item())
        D_losses.append(loss_D.item())

    print(
        f"Epoch [{epoch}/{n_epochs}] Loss_D: {loss_D.item()} Loss_G: {loss_G.item()}")

# 绘制 GP 损失曲线
plt.figure(figsize=(10, 5))
plt.plot(G_losses, label="Generator Loss")
plt.plot(D_losses, label="Discriminator Loss")
plt.legend()
plt.title("WGAN-GP Training Loss")
plt.xlabel("Iterations")
plt.ylabel("Loss")
plt.show()

# 绘制 DN 损失曲线
plt.figure(figsize=(10, 5))
plt.plot(D_losses, label="Discriminator Loss")
plt.legend()
plt.title("WGAN-GP Training Loss")
plt.xlabel("Iterations")
plt.ylabel("Loss")
plt.show()
