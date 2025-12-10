import tensorflow as tf
from tensorflow.keras import layers, Model

def conv_block(inputs, num_filters):
    """一个标准的卷积块"""
    x = layers.Conv2D(num_filters, 3, padding='same', activation='relu')(inputs)
    x = layers.Conv2D(num_filters, 3, padding='same', activation='relu')(x)
    return x

def encoder_block(inputs, num_filters):
    """编码器块 (卷积 + 下采样)"""
    cb = conv_block(inputs, num_filters)
    pool = layers.MaxPooling2D(2)(cb)
    return cb, pool

def decoder_block(inputs, skip_features, num_filters):
    """解码器块 (上采样 + 跳跃连接 + 卷积)"""
    x = layers.Conv2DTranspose(num_filters, 2, strides=2, padding='same')(inputs)
    x = layers.Concatenate()([x, skip_features])
    x = conv_block(x, num_filters)
    return x

def build_unet(input_shape=(32, 32, 1)):
    """构建一个 32x32 -> 32x32 的 U-Net 模型"""
    inputs = layers.Input(shape=input_shape)

    # 编码器 (下采样路径)
    s1, p1 = encoder_block(inputs, 16)
    s2, p2 = encoder_block(p1, 32)
    s3, p3 = encoder_block(p2, 64)
    
    # 瓶颈
    b1 = conv_block(p3, 128)
    
    # 解码器 (上采样路径)
    d1 = decoder_block(b1, s3, 64)
    d2 = decoder_block(d1, s2, 32)
    d3 = decoder_block(d2, s1, 16)
    
    # 输出层
    # 使用 'linear' 激活函数，因为我们预测的是强度值，而不是0-1的概率
    outputs = layers.Conv2D(1, 1, padding='same', activation='linear')(d3)

    model = Model(inputs, outputs, name='Sim-to-Real-UNet')
    return model

if __name__ == "__main__":
    # 快速测试模型结构
    model = build_unet()
    model.summary()

