#include <opencv2/opencv.hpp>
#include "minizbar.h"
#include <string>
#include <iostream>

using namespace cv;
namespace zbar = minizbar;


void zbarDecode(const cv::Mat& gray, std::string& code, cv::Mat& qr_coordinate) noexcept
{
    zbar::ImageScanner scanner; // 实例化扫描器
    scanner.set_config(zbar::ZBAR_NONE, zbar::ZBAR_CFG_ENABLE, 1);
    const uchar* raw = static_cast<uchar*>(gray.data);
    zbar::Image imageZbar(gray.cols, gray.rows, "Y800", raw, gray.cols * gray.rows);

    try {
        // 扫码
        scanner.scan(imageZbar);
        zbar::Image::SymbolIterator symbol = imageZbar.symbol_begin();
        if (imageZbar.symbol_begin() == imageZbar.symbol_end()) {
            std::cout << "zbar scan failed" << std::endl;
        }
        else {
            for (; symbol != imageZbar.symbol_end(); ++symbol) {
                code = symbol->get_data();
                qr_coordinate.at<cv::Point2f>(0) = cv::Point2f{ static_cast<float>(symbol->get_location_x(0)),static_cast<float>(symbol->get_location_y(0)) }; // 左下角
                qr_coordinate.at<cv::Point2f>(1) = cv::Point2f{ static_cast<float>(symbol->get_location_x(1)),static_cast<float>(symbol->get_location_y(1)) }; // 左上角
                qr_coordinate.at<cv::Point2f>(2) = cv::Point2f{ static_cast<float>(symbol->get_location_x(2)),static_cast<float>(symbol->get_location_y(2)) }; // 右上角
                qr_coordinate.at<cv::Point2f>(3) = cv::Point2f{ static_cast<float>(symbol->get_location_x(3)),static_cast<float>(symbol->get_location_y(3)) }; // 右下角
                break;
            }
        }
    }
    catch (const std::exception& e) {
        std::cout << "zbar decode error: " << e.what() << std::endl;
    }
    imageZbar.set_data(nullptr, 0); // 清理缓存
}

int main() {
    Mat gray = imread("D:/Donwloads/1.png", IMREAD_GRAYSCALE);
    cv::Mat qrs(4, 1, CV_32FC2);
    std::string s;
    zbarDecode(gray,s,qrs);

    std::cout << "s: " << s << ", qrs:" << qrs << std::endl;

	return 0;
}